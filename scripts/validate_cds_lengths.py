#!/usr/bin/env python3
"""
validate_cds_lengths.py

Check CDS length and basic sequence validity before Ka/Ks analysis.
"""

import argparse
from pathlib import Path


def read_fasta(path):
    seqs = {}
    current = None
    buf = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if current:
                    seqs[current] = "".join(buf).upper()
                current = s[1:].split()[0]
                buf = []
            else:
                buf.append(s)
        if current:
            seqs[current] = "".join(buf).upper()
    return seqs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cds", required=True)
    args = parser.parse_args()

    seqs = read_fasta(args.cds)
    print("seq_id\tlength\tmultiple_of_3\tstarts_with_ATG\tends_with_stop\tinvalid_bases")

    for seq_id, seq in seqs.items():
        invalid = sorted(set(seq) - set("ATCGN"))
        starts = seq.startswith("ATG")
        ends = seq[-3:] in {"TAA", "TAG", "TGA"} if len(seq) >= 3 else False
        print(
            f"{seq_id}\t{len(seq)}\t{len(seq) % 3 == 0}\t{starts}\t{ends}\t{''.join(invalid) if invalid else '-'}"
        )


if __name__ == "__main__":
    main()
