#!/usr/bin/env python3
# PURPOSE: dominant-chromosome PUR analysis (strict definition)
# INPUT:   nucmer coords
# OUTPUT:  PUR_strict results
# TOOL:    custom



import argparse
import os
import re
import sys


REF_CHR = {
    "NW_023276806.1": "chr1", "NW_023276807.1": "chr1",
    "NC_048595.1": "chr2", "NC_048596.1": "chr3", "NC_048597.1": "chr4",
    "NC_048598.1": "chr5", "NC_048599.1": "chr6", "NC_048600.1": "chr7",
    "NC_048601.1": "chr8", "NC_048602.1": "chr9", "NC_048603.1": "chr10",
    "NC_048604.1": "chrX",
}


def parse_coords(path, minlen):
    rows = []
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 13:
                continue
            try:
                len1 = int(c[4])
                s2, e2 = int(c[2]), int(c[3])
            except ValueError:
                continue
            if len1 < minlen:
                continue
            ref, qry = c[11], c[12]
            a, b = min(s2, e2), max(s2, e2)
            rows.append((qry, REF_CHR.get(ref, "unplaced"), a, b, len1))
    return rows


def merge_len(inter):
    inter = sorted(inter)
    if not inter:
        return 0
    tot = 0
    s, e = inter[0]
    for a, b in inter[1:]:
        if a <= e:
            e = max(e, b)
        else:
            tot += e - s + 1
            s, e = a, b
    tot += e - s + 1
    return tot


def complement(qlen, intervals):
    intervals = sorted(intervals)
    res = []
    cur = 1
    for a, b in intervals:
        if a > cur:
            res.append((cur, a - 1))
        cur = max(cur, b + 1)
    if cur <= qlen:
        res.append((cur, qlen))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coords1to1", required=True)
    ap.add_argument("--qlen", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--domminlen", type=int, default=100000)
    ap.add_argument("--minlen", type=int, default=2000)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    qlen = {}
    for line in open(args.qlen):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            qlen[c[0]] = int(c[1])

    rows = parse_coords(args.coords1to1, 1)
    print("alignments (LEN1>=1):", len(rows))

    dom_agg = {}
    for qry, rc, a, b, len1 in rows:
        if len1 >= args.domminlen:
            dom_agg.setdefault(qry, {}).setdefault(rc, 0)
            dom_agg[qry][rc] += len1

    dominant = {}
    for qry, d in dom_agg.items():
        dominant[qry] = max(d.items(), key=lambda kv: kv[1])[0]

    cov_inter = {}
    for qry, rc, a, b, len1 in rows:
        if len1 < args.minlen:
            continue
        if qry in dominant and rc == dominant[qry]:
            cov_inter.setdefault(qry, []).append((a, b))

    pur_bed = []
    per = []
    for s, L in sorted(qlen.items()):
        dom = dominant.get(s, "unplaced")
        cov = merge_len(cov_inter.get(s, []))
        pur = max(L - cov, 0)
        per.append((s, L, dom, cov, pur, pur / L * 100 if L else 0))
        for a, b in complement(L, cov_inter.get(s, [])):
            pur_bed.append((s, a - 1, b))

    with open(os.path.join(args.outdir, "pur_dominant.region"), "w") as f:
        for s, a, b in pur_bed:
            f.write(f"{s}\t{a}\t{b}\n")

    with open(os.path.join(args.outdir, "perChr_PUR_dominant.tsv"), "w") as f:
        f.write("scaffold\tlength\tdominant_ref\tcovered_bp\tPUR_bp\tPUR_pct\n")
        for s, L, dom, cov, pur, p in per:
            f.write(f"{s}\t{L}\t{dom}\t{cov}\t{pur}\t{p:.2f}\n")

    tot = sum(qlen.values())
    tot_cov = sum(merge_len(cov_inter.get(s, [])) for s in qlen)
    tot_pur = tot - tot_cov

    with open(os.path.join(args.outdir, "summary_pur_dominant.txt"), "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"query_scaffolds\t{len(qlen)}\n")
        f.write(f"query_total_bp\t{tot}\n")
        f.write(f"query_covered_to_dominant_bp\t{tot_cov}\n")
        f.write(f"query_covered_to_dominant_pct\t{tot_cov/tot*100:.2f}\n")
        f.write(f"PUR_dominant_bp\t{tot_pur}\n")
        f.write(f"PUR_dominant_pct\t{tot_pur/tot*100:.2f}\n")
        f.write(f"PUR_dominant_ge1kb_bp\t{sum(b-a for _,a,b in pur_bed if b-a>=1000)}\n")
        f.write(f"PUR_dominant_ge1kb_count\t{sum(1 for _,a,b in pur_bed if b-a>=1000)}\n")
        no_dom = sum(1 for s in qlen if dominant.get(s, "unplaced") == "unplaced")
        f.write(f"scaffolds_dominant_unplaced\t{no_dom}\n")

    print(f"PUR_dominant = {tot_pur:,} bp ({tot_pur/tot*100:.2f}% of {tot:,})")
    print(f"covered to dominant chr = {tot_cov:,} bp ({tot_cov/tot*100:.2f}%)")
    print("outputs ->", args.outdir)


if __name__ == "__main__":
    main()
