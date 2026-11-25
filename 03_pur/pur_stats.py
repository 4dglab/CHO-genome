#!/usr/bin/env python3
# PURPOSE: PUR statistics (covered / uncovered intervals)
# INPUT:   nucmer coords
# OUTPUT:  PUR intervals (BED)
# TOOL:    custom



import argparse
import os
import re
import sys


def parse_coords(path):
    rows = []
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 13:
                continue
            try:
                s1, e1, s2, e2 = (int(x) for x in c[:4])
                len1, len2 = (int(x) for x in c[4:6])
                lenr, lenq = (int(x) for x in c[7:9])
            except ValueError:
                continue
            ref, qry = c[11], c[12]
            strand = "-" if s2 > e2 else "+"
            rows.append({
                "ref": ref, "qry": qry,
                "rs": min(s1, e1), "re": max(s1, e1),
                "qs": min(s2, e2), "qe": max(s2, e2),
                "strand": strand, "len1": len1, "len2": len2,
            })
    return rows


def merge_len(intervals):
    
    if not intervals:
        return 0
    intervals = sorted(intervals)
    tot = 0
    s, e = intervals[0]
    for a, b in intervals[1:]:
        if a <= e:
            e = max(e, b)
        else:
            tot += e - s + 1
            s, e = a, b
    tot += e - s + 1
    return tot


def clip(intervals, length):
    
    out = []
    for a, b in intervals:
        a = max(a, 1)
        b = min(b, length)
        if a <= b:
            out.append((a, b))
    return out


def complement(qlen, intervals):
    
    if not intervals:
        return [(1, qlen)]
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
    ap.add_argument("--coords", required=True)
    ap.add_argument("--qlen", required=True)
    ap.add_argument("--rlen", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--minlen", type=int, default=2000,
                    help="minimum alignment length (LEN1) to keep (default 2000, ~ delta-filter -l)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    qlen = {}
    for line in open(args.qlen):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            qlen[c[0]] = int(c[1])
    rlen = {}
    for line in open(args.rlen):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            rlen[c[0]] = int(c[1])

    rows = [r for r in parse_coords(args.coords) if r["len1"] >= args.minlen]
    print("alignments parsed:", len(rows), f"(LEN1 >= {args.minlen})")

    q_intervals = {}
    r_intervals = {}
    for r in rows:
        q_intervals.setdefault(r["qry"], []).append((r["qs"], r["qe"]))
        r_intervals.setdefault(r["ref"], []).append((r["rs"], r["re"]))

    q_cov = {s: merge_len(clip(v, qlen[s])) for s, v in q_intervals.items()}

    # per-scaffold PUR (BED 0-based)
    pur_bed = []
    per_chr = []
    for s, L in sorted(qlen.items()):
        cov = q_cov.get(s, 0)
        pur = max(L - cov, 0)
        per_chr.append((s, L, cov, pur, pur / L * 100 if L else 0))
        for a, b in complement(L, clip(q_intervals.get(s, []), L)):
            pur_bed.append((s, a - 1, b))

    with open(os.path.join(args.outdir, "pur_full.region"), "w") as f:
        for s, a, b in pur_bed:
            f.write(f"{s}\t{a}\t{b}\n")

    with open(os.path.join(args.outdir, "perChr_PUR.tsv"), "w") as f:
        f.write("scaffold\tlength\tcovered_bp\tPUR_bp\tPUR_pct\n")
        for s, L, cov, pur, p in per_chr:
            f.write(f"{s}\t{L}\t{cov}\t{pur}\t{p:.2f}\n")

    small = [(s, L, cov, pur, cov / L * 100 if L else 0)
             for s, L, cov, pur, p in per_chr
             if re.fullmatch(r"scaffold_(\d+)", s) and int(re.match(r"scaffold_(\d+)", s).group(1)) > 11]

    def cat(pct):
        if pct >= 50:
            return "anchored_to_ref"
        if pct >= 10:
            return "partial"
        return "novel_or_unplaced"

    with open(os.path.join(args.outdir, "small_scaffolds_classify.tsv"), "w") as f:
        f.write("scaffold\tlength\taligned_bp\taligned_pct\tcategory\n")
        for s, L, cov, pur, p in small:
            f.write(f"{s}\t{L}\t{cov}\t{p:.2f}\t{cat(p)}\n")

    r_cov = {s: merge_len(clip(v, rlen[s])) for s, v in r_intervals.items()}
    ref_total = sum(rlen.values())
    ref_cov_total = sum(r_cov.get(s, 0) for s in rlen)
    ref_uncovered = max(ref_total - ref_cov_total, 0)

    q_total = sum(qlen.values())
    q_cov_total = sum(q_cov.get(s, 0) for s in qlen)
    q_pur = max(q_total - q_cov_total, 0)

    small_anch = sum(1 for _, _, _, _, p in small if cat(p) == "anchored_to_ref")
    small_novel = sum(1 for _, _, _, _, p in small if cat(p) == "novel_or_unplaced")
    small_part = len(small) - small_anch - small_novel

    with open(os.path.join(args.outdir, "summary_pur_full.txt"), "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"query_scaffolds\t{len(qlen)}\n")
        f.write(f"query_total_bp\t{q_total}\n")
        f.write(f"query_covered_bp\t{q_cov_total}\n")
        f.write(f"query_covered_pct\t{q_cov_total/q_total*100:.2f}\n")
        f.write(f"PUR_full_bp\t{q_pur}\n")
        f.write(f"PUR_full_pct\t{q_pur/q_total*100:.2f}\n")
        f.write(f"PUR_full_ge1kb_bp\t{sum(b-a for _,a,b in pur_bed if b-a>=1000)}\n")
        f.write(f"PUR_full_ge1kb_count\t{sum(1 for _,a,b in pur_bed if b-a>=1000)}\n")
        f.write(f"ref_sequences\t{len(rlen)}\n")
        f.write(f"ref_total_bp\t{ref_total}\n")
        f.write(f"ref_covered_bp\t{ref_cov_total}\n")
        f.write(f"ref_covered_pct\t{ref_cov_total/ref_total*100:.2f}\n")
        f.write(f"ref_uncovered_bp\t{ref_uncovered}\n")
        f.write(f"ref_uncovered_pct\t{ref_uncovered/ref_total*100:.2f}\n")
        f.write(f"small_scaffolds_total\t{len(small)}\n")
        f.write(f"small_anchored_to_ref\t{small_anch}\n")
        f.write(f"small_partial\t{small_part}\n")
        f.write(f"small_novel_or_unplaced\t{small_novel}\n")

    print("=== summary ===")
    print(f"query covered {q_cov_total/q_total*100:.2f}%  -> PUR_full {q_pur} bp ({q_pur/q_total*100:.2f}%)")
    print(f"ref covered {ref_cov_total/ref_total*100:.2f}%  -> ref uncovered {ref_uncovered} bp")
    print(f"small scaffolds: anchored {small_anch}, partial {small_part}, novel {small_novel}")
    print("outputs ->", args.outdir)


if __name__ == "__main__":
    main()
