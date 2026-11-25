#!/usr/bin/env python3
# PURPOSE: Final GFF + filtering pipeline report
# INPUT:   GFF
# OUTPUT:  final GFF + filtering_pipeline.tsv
# TOOL:    custom


"""Filtering pipeline report + final GFF construction (B4 step).

Reads the evidence-annotated GFF and the two drop lists (pseudogene + redundant),
produces:
  final/genomic.gff                retained genes (all isoforms), Evidence_level kept
  final/kept_genes.tsv             retained gene ids
  final/filtering_pipeline.tsv     stage | genes | mRNAs | transcripts | dropped | note
  final/filtering_pipeline.png     bar chart (waterfall style)

Usage:
  filtering_report.py --gff <evidence/genomic.evidence.gff> \
      --pseudo-drop <pseudogene/pseudogene_dropped.tsv> \
      --redundant-drop <redundancy/dropped_redundant.tsv> \
      --outdir <final/>
"""
import argparse
import os
import re
from collections import Counter


def load_dropped(path):
    ids = set()
    if path and os.path.exists(path):
        for line in open(path):
            if line.startswith("gene_id"):
                continue
            f = line.rstrip("\n").split("\t")
            if f:
                ids.add(f[0])
    return ids


def count_features(gff_path):
    genes = mRNAs = transcripts = 0
    for line in open(gff_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        if f[2] == "gene":
            genes += 1
        elif f[2] == "mRNA":
            mRNAs += 1
        elif f[2] == "transcript":
            transcripts += 1
    return genes, mRNAs, transcripts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--pseudo-drop", required=True)
    ap.add_argument("--redundant-drop", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pseudo_drop = load_dropped(args.pseudo_drop)
    redun_drop = load_dropped(args.redundant_drop)
    all_drop = pseudo_drop | redun_drop

    # Pre-pass: map every transcript (and thus its children) to its root gene,
    # so that dropping a gene also removes its mRNA/exon/CDS children.
    # (Fix: previously children were compared against gene-level drop lists
    #  via their DIRECT parent (transcript id) and therefore survived as
    #  orphaned features.)
    tx2gene = {}
    with open(args.gff) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(re.findall(r"([^=;]+)=([^;]+)", f[8]))
            fid = attrs.get("ID")
            par = attrs.get("Parent")
            if f[2] in ("mRNA", "transcript", "lnc_RNA") and fid and par:
                tx2gene[fid] = par.split(",")[0]

    def root_gene(ftype, attrs):
        if ftype == "gene":
            return attrs.get("ID")
        pid = attrs.get("Parent", "").split(",")[0]
        if ftype in ("mRNA", "transcript", "lnc_RNA"):
            return pid
        return tx2gene.get(pid)

    # final GFF: keep gene + children unless ROOT gene dropped
    final_gff = os.path.join(args.outdir, "genomic.gff")
    kept_ids = set()
    with open(args.gff) as fin, open(final_gff, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                fout.write(line)
                continue
            a = dict(re.findall(r"([^=;]+)=([^;]+)", f[8]))
            root = root_gene(f[2], a)
            if root is not None and root in all_drop:
                continue
            if f[2] == "gene":
                kept_ids.add(root)
            fout.write(line)

    with open(os.path.join(args.outdir, "kept_genes.tsv"), "w") as out:
        for g in sorted(kept_ids):
            out.write(g + "\n")

    # ---- pipeline counts ----
    initial_genes, initial_mrna, _ = count_features(args.gff)
    after_pseudo = os.path.join(args.outdir, "stage_after_pseudo.gff")
    pseudo_kept_ids = kept_ids | redun_drop  # genes surviving pseudo step
    # retained if not in pseudo_drop
    with open(args.gff) as fin, open(after_pseudo, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                fout.write(line); continue
            a = dict(re.findall(r"([^=;]+)=([^;]+)", f[8]))
            root = root_gene(f[2], a)
            if root is not None and root in pseudo_drop:
                continue
            fout.write(line)
    p_g, p_m, p_t = count_features(after_pseudo)
    f_g, f_m, f_t = count_features(final_gff)

    rows = [
        ("initial (after AGAT UTR filter)", initial_genes, initial_mrna, 0,
         "", ""),
        ("after pseudogene filtering", p_g, p_m, initial_genes - p_g,
         "pseudogene/truncated", "%d genes dropped" % (initial_genes - p_g)),
        ("after redundancy removal", f_g, f_m, p_g - f_g,
         "redundant/nested/CDS-overlap", "%d genes dropped" % (p_g - f_g)),
        ("final retained", f_g, f_m, 0, "", ""),
    ]
    with open(os.path.join(args.outdir, "filtering_pipeline.tsv"), "w") as out:
        out.write("stage\tgenes\tmRNAs\tdropped\tcategory\tnote\n")
        for r in rows:
            out.write("%s\t%d\t%d\t%d\t%s\t%s\n" % r)

    # ---- bar chart (matplotlib, base env) ----
    if not args.no_figure:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available; skipping figure")
            return
        stages = [r[0] for r in rows]
        counts = [r[1] for r in rows]
        colors = ["#4C72B0", "#DD8452", "#DD8452", "#55A868"]
        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(stages, counts, color=colors, edgecolor="black", linewidth=0.8)
        ax.set_ylabel("Number of genes")
        ax.set_ylim(0, max(counts) * 1.15)
        for b, c in zip(bars, counts):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 5,
                    str(c), ha="center", va="bottom", fontsize=10)
        ax.set_title("Gene filtering pipeline")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "filtering_pipeline.png"), dpi=200)
        print("figure written")

    print("FILTERING_REPORT_DONE")
    print("  initial genes: %d  mRNAs: %d" % (initial_genes, initial_mrna))
    print("  after pseudogene: %d genes" % p_g)
    print("  final retained: %d genes  %d mRNAs  %d augustus transcripts"
          % (f_g, f_m, f_t))
    print("  pseudogene drops: %d  redundant drops: %d" %
          (len(pseudo_drop), len(redun_drop)))


if __name__ == "__main__":
    main()
