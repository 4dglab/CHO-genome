#!/usr/bin/env python3
# PURPOSE: BLAST conservation statistics
# INPUT:   blast -outfmt 6
# OUTPUT:  pur_vs_species_stats.txt
# TOOL:    custom



import argparse
import os


def parse_tsv(path):
    hits = {}
    if not path or not os.path.exists(path):
        return hits
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) < 12:
            continue
        try:
            qseqid, qlen = c[0], int(c[1])
            pident, length = float(c[4]), int(c[5])
        except ValueError:
            continue
        if length >= 100 and pident >= 80:
            if qseqid not in hits or hits[qseqid][1] < length:
                hits[qseqid] = (qlen, length, pident)
    return hits


def fasta_len(path):
    n = b = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                n += 1
            else:
                b += len(line.strip())
    return n, b


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mouse", default=None)
    ap.add_argument("--rat", default=None)
    ap.add_argument("--purfa", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    n_pur, b_pur = fasta_len(args.purfa)
    print(f"PUR >=1kb: {n_pur} seqs, {b_pur} bp")

    with open(args.out, "w") as f:
        f.write("species\tPUR_seqs_with_hit\tPUR_seqs_pct\tPUR_bp_with_hit\tPUR_bp_pct\tmean_pident\tmedian_len\n")
        for label, p in (("mouse", args.mouse), ("rat", args.rat)):
            hits = parse_tsv(p)
            if not hits:
                f.write(f"{label}\t0\t0\t0\t0\t-\t-\n")
                print(f"{label}: no hits (file missing or filtered out)")
                continue
            seqs = len(hits)
            bp = sum(v[0] for v in hits.values())
            pids = sorted(v[2] for v in hits.values())
            lens = sorted(v[1] for v in hits.values())
            mean_pid = sum(pids) / len(pids)
            med_len = lens[len(lens) // 2]
            f.write(f"{label}\t{seqs}\t{seqs/n_pur*100:.2f}\t{bp}\t{bp/b_pur*100:.2f}\t{mean_pid:.2f}\t{med_len}\n")
            print(f"{label}: {seqs} seqs ({seqs/n_pur*100:.2f}%), {bp} bp ({bp/b_pur*100:.2f}%), "
                  f"mean pident {mean_pid:.2f}, median len {med_len}")

    print("output ->", args.out)


if __name__ == "__main__":
    main()
