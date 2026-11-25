#!/usr/bin/env python3
# PURPOSE: Soft-mask genome by repeat annotation
# INPUT:   genome.fa + repeat.gff3
# OUTPUT:  genome.softmasked.fa
# TOOL:    custom


"""Soft-mask a genome FASTA using intervals from a GFF3 file.

Repeats (from merge_all.repeat.gff3) are lowercased so Augustus --softmasking=1
ignores them during ab initio prediction.

Usage:
    softmask_genome.py genome.fa repeats.gff3 out.softmasked.fa
"""

import sys
from collections import defaultdict


def parse_gff3_intervals(gff_path):
    intervals = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            seqid, start, end = cols[0], int(cols[3]), int(cols[4])
            if start > end:
                start, end = end, start
            intervals[seqid].append((start, end))
    return intervals


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs.sort()
    merged = []
    cur_s, cur_e = ivs[0]
    for s, e in ivs[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    genome_fa, gff, out_fa = sys.argv[1], sys.argv[2], sys.argv[3]

    intervals = {s: merge_intervals(iv) for s, iv in parse_gff3_intervals(gff).items()}

    n_bp = 0
    with open(genome_fa) as fh, open(out_fa, "w") as out:
        seqid = None
        seq_parts = []
        for line in fh:
            if line.startswith(">"):
                if seqid is not None:
                    n_bp += flush_seq(seqid, seq_parts, intervals.get(seqid, []), out)
                seqid = line.split()[0][1:]
                seq_parts = []
                out.write(line)
                continue
            if line.strip() == "":
                continue
            seq_parts.append(line.rstrip("\n"))
        if seqid is not None:
            n_bp += flush_seq(seqid, seq_parts, intervals.get(seqid, []), out)

    print(f"soft-masked: {n_bp} bp total lowercase across genome", file=sys.stderr)


def flush_seq(seqid, seq_parts, iv_list, out):
    seq = bytearray("".join(seq_parts).encode("ascii"))
    iv_list = merge_intervals(iv_list)
    masked = 0
    for s, e in iv_list:
        s = max(1, s)
        e = min(len(seq), e)
        if s > e:
            continue
        for i in range(s - 1, e):
            c = seq[i]
            if 65 <= c <= 90:
                seq[i] = c + 32
                masked += 1
    for i in range(0, len(seq), 80):
        out.write(seq[i : i + 80].decode("ascii") + "\n")
    return masked


if __name__ == "__main__":
    main()
