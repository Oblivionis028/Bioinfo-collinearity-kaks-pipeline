#!/usr/bin/env python3
"""
prepare_mcscanx_inputs.py

Generate MCScanX-compatible input files from a GFF/GFF3 annotation file and
a protein FASTA file.

Outputs:
  <prefix>.clean.protein.fasta
  <prefix>.gff

The MCScanX 4-column GFF format is:
  chromosome/scaffold    gene_id/protein_id    start    end

No real data is included in this script.
"""

import argparse
import re
from pathlib import Path


def clean_fasta_header(protein_path: Path, out_path: Path):
    ids = []
    with protein_path.open("r", encoding="utf-8", errors="replace") as fin, out_path.open("w", encoding="utf-8", newline="\n") as fout:
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                ids.append(seq_id)
                fout.write(f">{seq_id}\n")
            else:
                fout.write(line.rstrip("\n") + "\n")
    return ids


def get_attr(attrs: str, key: str):
    pattern = rf"(?:^|;){re.escape(key)}=([^;]+)"
    m = re.search(pattern, attrs)
    return m.group(1) if m else None


def parse_gff_cds_by_protein(gff_path: Path, feature_type: str, protein_attr: str):
    coords = {}

    with gff_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _source, feature, start, end, _score, _strand, _phase, attrs = parts
            if feature != feature_type:
                continue

            gene_id = get_attr(attrs, protein_attr)
            if not gene_id:
                continue

            try:
                s, e = int(start), int(end)
            except ValueError:
                continue

            if gene_id not in coords:
                coords[gene_id] = [chrom, s, e]
            else:
                coords[gene_id][1] = min(coords[gene_id][1], s)
                coords[gene_id][2] = max(coords[gene_id][2], e)

    return coords


def write_mcscanx_gff(coords, fasta_ids, out_path: Path):
    written = 0
    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        for gene_id in fasta_ids:
            if gene_id in coords:
                chrom, start, end = coords[gene_id]
                out.write(f"{chrom}\t{gene_id}\t{start}\t{end}\n")
                written += 1
    return written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, help="Input GFF/GFF3 file")
    parser.add_argument("--protein", required=True, help="Input protein FASTA file")
    parser.add_argument("--prefix", default="sample", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--feature_type", default="CDS", help="GFF feature type used to define coordinates")
    parser.add_argument("--protein_attr", default="Protein_Accession", help="GFF attribute key matching protein FASTA IDs")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gff_path = Path(args.gff)
    protein_path = Path(args.protein)

    clean_protein = outdir / f"{args.prefix}.clean.protein.fasta"
    mcscanx_gff = outdir / f"{args.prefix}.gff"

    fasta_ids = clean_fasta_header(protein_path, clean_protein)
    coords = parse_gff_cds_by_protein(gff_path, args.feature_type, args.protein_attr)
    written = write_mcscanx_gff(coords, fasta_ids, mcscanx_gff)

    print("Done.")
    print(f"Protein records in FASTA: {len(fasta_ids)}")
    print(f"Protein IDs with coordinates in GFF: {len(coords)}")
    print(f"Rows written to {mcscanx_gff}: {written}")
    print()
    print("Generated files:")
    print(f"  {clean_protein}")
    print(f"  {mcscanx_gff}")
    print()
    print("MCScanX GFF column order:")
    print("  chromosome/scaffold    gene_id/protein_id    start    end")

    if written == 0:
        print()
        print("ERROR: No rows were written.")
        print(f"Check whether GFF feature '{args.feature_type}' contains attribute '{args.protein_attr}=...'")
    elif written != len(fasta_ids):
        print()
        print("WARNING: Some protein FASTA IDs were not found in the GFF-derived coordinates.")


if __name__ == "__main__":
    main()
