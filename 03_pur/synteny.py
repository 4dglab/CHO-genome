#!/usr/bin/env python3
# PURPOSE: synteny block analysis
# INPUT:   nucmer coords
# OUTPUT:  synteny results
# TOOL:    custom



import argparse
import os
import re
import statistics


def ref_chr_map(ref_fa):
    m = {}
    for line in open(ref_fa):
        if not line.startswith(">"):
            continue
        hdr = line[1:].strip()
        name = hdr.split()[0]
        if "chr1_0" in hdr or "chr1_1" in hdr:
            m[name] = "chr1"
        elif "chromosome X" in hdr:
            m[name] = "chrX"
        else:
            mo = re.search(r"chromosome (\d+)", hdr)
            m[name] = "chr" + mo.group(1) if mo else "unplaced"
    return m


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
            except ValueError:
                continue
            ref, qry = c[11], c[12]
            rows.append({
                "ref": ref, "qry": qry,
                "rs": min(s1, e1), "re": max(s1, e1),
                "qs": min(s2, e2), "qe": max(s2, e2),
                "strand": "-" if s2 > e2 else "+",
                "len1": len1, "len2": len2,
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coords", required=True)
    ap.add_argument("--ref-fa", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--minlen", type=int, default=2000,
                    help="minimum alignment length (LEN1) to keep (default 2000)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rcmap = ref_chr_map(args.ref_fa)
    rows = [r for r in parse_coords(args.coords) if r["len1"] >= args.minlen]
    print("alignments parsed:", len(rows), f"(LEN1 >= {args.minlen})")

    q_chrom = {f"scaffold_{i}": ("chr%d" % i if i <= 10 else "chrX")
               for i in range(1, 12)}

    agg = {}  # qry -> {refchr: bp}
    for r in rows:
        rc = rcmap.get(r["ref"], "unplaced")
        agg.setdefault(r["qry"], {}).setdefault(rc, 0)
        agg[r["qry"]][rc] += r["len1"]

    # dominant assignment
    dom_rows = []
    for qry, d in agg.items():
        ranked = sorted(d.items(), key=lambda kv: kv[1], reverse=True)
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else (None, 0)
        total = sum(d.values())
        frac = best[1] / total * 100 if total else 0
        dom_rows.append((qry, best[0], best[1], total, frac,
                         (second[0] + ":" + str(second[1])) if second[0] else "-"))

    with open(os.path.join(args.outdir, "scaffold_dominant.tsv"), "w") as f:
        f.write("scaffold\tdominant_ref_chr\tdominant_bp\ttotal_aligned_bp\tdominant_pct\tsecond\n")
        for qry, best, bb, tot, frac, sec in sorted(dom_rows, key=lambda x: x[3], reverse=True):
            f.write(f"{qry}\t{best}\t{bb}\t{tot}\t{frac:.2f}\t{sec}\n")

    blocks = []
    for r in rows:
        if r["qry"] in q_chrom:
            rc = rcmap.get(r["ref"], "unplaced")
            dom = dict((q, c) for q, c, *_ in dom_rows).get(r["qry"])
            if rc == dom:
                blocks.append(r)

    blocks.sort(key=lambda r: r["len1"], reverse=True)
    lens = [b["len1"] for b in blocks]
    inv = [b for b in blocks if b["strand"] == "-"]
    inv_bp = sum(b["len1"] for b in inv)

    with open(os.path.join(args.outdir, "blocks.strand.tsv"), "w") as f:
        f.write("refchr\tqrychr\trstart\trend\tqstart\tqend\tlen\tstrand\n")
        for b in blocks:
            rc = rcmap.get(b["ref"])
            f.write(f"{rc}\t{q_chrom[b['qry']]}\t{b['rs']}\t{b['re']}\t{b['qs']}\t{b['qe']}\t{b['len1']}\t{b['strand']}\n")

    inter = []
    for qry, best, bb, tot, frac, sec in dom_rows:
        if qry in q_chrom:
            for rc, bp in agg[qry].items():
                if rc != best:
                    inter.append((qry, q_chrom[qry], best, rc, bp))

    with open(os.path.join(args.outdir, "summary_synteny.txt"), "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"blocks_total\t{len(lens)}\n")
        f.write(f"blocks_median_len\t{statistics.median(lens) if lens else 0}\n")
        f.write(f"blocks_mean_len\t{statistics.mean(lens) if lens else 0:.0f}\n")
        f.write(f"blocks_N50_len\t{stats_n50(lens)}\n")
        f.write(f"blocks_total_bp\t{sum(lens)}\n")
        f.write(f"inversion_count\t{len(inv)}\n")
        f.write(f"inversion_pct_count\t{len(inv)/len(lens)*100 if lens else 0:.2f}\n")
        f.write(f"inversion_bp\t{inv_bp}\n")
        f.write(f"inversion_pct_bp\t{inv_bp/sum(lens)*100 if lens else 0:.2f}\n")
        f.write(f"interchromosomal_signals\t{len(inter)}\n")

    print("=== synteny summary ===")
    print(f"blocks {len(lens)}, N50 {stats_n50(lens)}, median {statistics.median(lens) if lens else 0}")
    print(f"inversions {len(inv)} ({len(inv)/len(lens)*100 if lens else 0:.2f}% by count), {inv_bp} bp")
    print("--- chromosome correspondence ---")
    for qry, best, bb, tot, frac, sec in sorted(dom_rows, key=lambda x: x[3], reverse=True):
        if qry in q_chrom:
            print(f"  {qry:12s} -> {best:5s} {frac:6.2f}%  (aligned {bb} bp)")
    print("--- inter-chromosomal signals ---")
    for qry, qc, best, rc, bp in inter:
        print(f"  {qry} ({qc}) has {bp} bp aligned to {rc} (dominant {best})")
    print("outputs ->", args.outdir)


def stats_n50(lens):
    if not lens:
        return 0
    half = sum(lens) / 2
    cum = 0
    for L in sorted(lens, reverse=True):
        cum += L
        if cum >= half:
            return L
    return 0


if __name__ == "__main__":
    main()
