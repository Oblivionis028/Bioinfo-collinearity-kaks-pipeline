#!/usr/bin/env python3
"""
filter_collinearity_targets.py

Filter target-related collinear gene pairs from an MCScanX .collinearity file.

Inputs:
  --collinearity  MCScanX .collinearity file
  --gff           MCScanX 4-column GFF: chromosome/scaffold gene_id start end
  --targets       one target ID per line

Outputs:
  <out_prefix>_pairs.tsv
  <out_prefix>_summary.tsv
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict


def read_targets(path: Path):
    targets = set()
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                targets.add(s.split()[0])
    return targets


def read_gff4(path: Path):
    loc = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 4:
                continue
            chrom, gene, start, end = parts[:4]
            try:
                loc[gene] = (chrom, int(start), int(end))
            except ValueError:
                continue
    return loc


def parse_collinearity(path: Path):
    pairs = []
    block_id = ""
    block_header = ""
    token_re = re.compile(r"[A-Za-z0-9_.:-]+")

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("## Alignment"):
                block_header = s
                m = re.search(r"Alignment\s+(\d+)", s)
                block_id = m.group(1) if m else ""
                continue

            if s.startswith("#"):
                continue

            tokens = token_re.findall(s)
            # MCScanX pair lines normally contain exactly two gene IDs after a position field.
            # To avoid assuming a specific ID prefix, take tokens that are not purely numeric
            # and are not obvious e-values/scores.
            gene_like = []
            for t in tokens:
                if re.fullmatch(r"\d+", t):
                    continue
                if re.fullmatch(r"\d+e[-+]?\d+", t, re.IGNORECASE):
                    continue
                if t in {"plus", "minus"}:
                    continue
                gene_like.append(t)

            # Robust fallback: MCScanX lines often look like:
            # 0- 0: geneA geneB evalue
            # After filtering, last two non-numeric IDs are usually geneA/geneB.
            if len(gene_like) >= 2:
                gene1, gene2 = gene_like[-2], gene_like[-1]
                fields = s.split()
                evalue = fields[-1] if fields else ""
                pairs.append({
                    "block_id": block_id,
                    "block_header": block_header,
                    "gene1": gene1,
                    "gene2": gene2,
                    "evalue_or_score": evalue,
                    "raw_line": s,
                })

    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--collinearity", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--out_prefix", default="target_collinearity")
    args = parser.parse_args()

    targets = read_targets(Path(args.targets))
    loc = read_gff4(Path(args.gff))
    all_pairs = parse_collinearity(Path(args.collinearity))

    target_rows = []
    summary = defaultdict(lambda: {
        "partners": set(),
        "blocks": set(),
        "same": 0,
        "different": 0,
        "total": 0,
    })

    for p in all_pairs:
        g1, g2 = p["gene1"], p["gene2"]
        hit_targets = []
        if g1 in targets:
            hit_targets.append((g1, g2))
        if g2 in targets:
            hit_targets.append((g2, g1))

        for target, partner in hit_targets:
            t_loc = loc.get(target, ("NA", "NA", "NA"))
            p_loc = loc.get(partner, ("NA", "NA", "NA"))
            relation = "same_scaffold_or_chromosome" if t_loc[0] == p_loc[0] and t_loc[0] != "NA" else "different_scaffold_or_chromosome"

            target_rows.append({
                "target_gene": target,
                "partner_gene": partner,
                "block_id": p["block_id"],
                "target_chr": t_loc[0],
                "target_start": t_loc[1],
                "target_end": t_loc[2],
                "partner_chr": p_loc[0],
                "partner_start": p_loc[1],
                "partner_end": p_loc[2],
                "relation": relation,
                "evalue_or_score": p["evalue_or_score"],
                "block_header": p["block_header"],
            })

            s = summary[target]
            s["partners"].add(partner)
            s["blocks"].add(p["block_id"])
            s["total"] += 1
            if relation == "same_scaffold_or_chromosome":
                s["same"] += 1
            else:
                s["different"] += 1

    pairs_out = Path(f"{args.out_prefix}_pairs.tsv")
    summary_out = Path(f"{args.out_prefix}_summary.tsv")

    pair_headers = [
        "target_gene", "partner_gene", "block_id",
        "target_chr", "target_start", "target_end",
        "partner_chr", "partner_start", "partner_end",
        "relation", "evalue_or_score", "block_header"
    ]
    with pairs_out.open("w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(pair_headers) + "\n")
        for r in sorted(target_rows, key=lambda x: (x["target_gene"], x["block_id"], x["partner_gene"])):
            out.write("\t".join(str(r[h]) for h in pair_headers) + "\n")

    summary_headers = [
        "target_gene", "found_in_collinearity",
        "num_partner_genes", "num_blocks", "num_pairs",
        "same_scaffold_or_chromosome_pairs",
        "different_scaffold_or_chromosome_pairs",
        "partner_genes"
    ]
    with summary_out.open("w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(summary_headers) + "\n")
        for target in sorted(targets):
            if target in summary:
                s = summary[target]
                row = [
                    target, "YES", len(s["partners"]), len(s["blocks"]), s["total"],
                    s["same"], s["different"], ";".join(sorted(s["partners"]))
                ]
            else:
                row = [target, "NO", 0, 0, 0, 0, 0, ""]
            out.write("\t".join(map(str, row)) + "\n")

    print("Done.")
    print(f"All collinear gene pairs parsed: {len(all_pairs)}")
    print(f"Target-related collinear rows: {len(target_rows)}")
    print(f"Targets found in collinearity: {sum(1 for t in targets if t in summary)} / {len(targets)}")
    print()
    print("Generated:")
    print(f"  {pairs_out}")
    print(f"  {summary_out}")


if __name__ == "__main__":
    main()
