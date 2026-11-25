#!/usr/bin/env python3
"""Encode OnTAD .tad files into a merged TSV with TAD position and nesting level.

Usage: 03_encode_tad.py [reslabel] [res_bp]
  defaults: reslabel=10kb, res_bp=10000
Output: ${HI_C_DIR}/TAD/combined/onTAD.<reslabel>.encoded.tsv
Columns: chrom start end tad_id TADlevel level TADmean TADscore
  level = clamped label (level1..level5); true TADlevel >=5 -> level5.
  The level-0 whole-chromosome root is excluded.

Env (set in config.sh):
  HI_C_DIR             output root of the Hi-C module (TAD/<reslabel>.onTAD/)
"""
import os
import sys
import glob

HI_C_DIR = os.environ.get("HI_C_DIR", os.path.join(os.getcwd(), "hi_c"))


def main(reslabel="10kb", res_bp=10000, outdir=None):
    base = os.path.join(HI_C_DIR, "TAD", f"{reslabel}.onTAD")
    outdir = outdir or os.path.join(HI_C_DIR, "TAD", "combined")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"onTAD.{reslabel}.encoded.tsv")
    rows = []
    for tad in sorted(glob.glob(os.path.join(base, "onTAD.K1.scaffold_*.tad"))):
        chrom = os.path.basename(tad).replace("onTAD.K1.", "").replace(".tad", "")
        idx = 0
        for line in open(tad):
            line = line.strip()
            if not line:
                continue
            p = line.split()
            if len(p) != 5:
                continue
            sb, e, tl, mean, score = int(p[0]), int(p[1]), int(p[2]), float(p[3]), float(p[4])
            if tl == 0:  # skip whole-chromosome root
                continue
            idx += 1
            start = sb * res_bp
            end = e * res_bp
            level = f"level{tl}" if tl <= 5 else "level5"
            tid = f"{chrom}_tad_{idx:06d}"
            rows.append([chrom, start, end, tid, tl, level, mean, score])
    with open(outpath, "w") as o:
        o.write("chrom\tstart\tend\ttad_id\tTADlevel\tlevel\tTADmean\tTADscore\n")
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")
    print(f"WRITTEN {len(rows)} TADs -> {outpath}")


if __name__ == "__main__":
    rl = sys.argv[1] if len(sys.argv) > 1 else "10kb"
    rb = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    od = sys.argv[3] if len(sys.argv) > 3 else None
    main(rl, rb, od)
