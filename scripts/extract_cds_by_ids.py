#!/usr/bin/env python3
"""
extract_cds_by_ids.py

Extract selected CDS sequences from a FASTA file.

Inputs:
  --cds  all CDS FASTA
  --ids  one sequence ID per line
  --out  output FASTA
"""

import argparse
from pathlib import Path


def read_ids(path: Path):
    ids = set()
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.add(s.split()[0])
    return ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cds", required=True)
    parser.add_argument("--ids", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    wanted = read_ids(Path(args.ids))
    found = set()

    out = open(args.out, "w", encoding="utf-8", newline="\n")
    write = False

    with open(args.cds, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                write = seq_id in wanted
                if write:
                    found.add(seq_id)
                    out.write(line)
            else:
                if write:
                    out.write(line)

    out.close()

    print(f"Found: {len(found)} / {len(wanted)}")
    missing = sorted(wanted - found)
    print(f"Missing: {missing}")
    print(f"Output: {args.out}")


if __name__ == "__main__":
    main()
