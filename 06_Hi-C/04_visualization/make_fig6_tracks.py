#!/usr/bin/env python3
# env: hic-viz-env
# make_fig6_tracks.py — generate Figure6-style tracks for given genes/loci,
# replicating the two notebooks:
#   Figure6_2genes_trackc_newref.ipynb    -> <GENE>_state_trackc.pdf
#   Figure6_markline_2genes_newref.ipynb  -> <GENE>_markline.pdf/.png
#
# Usage: python3 make_fig6_tracks.py
# Outputs are written under result/revision/bigwig/figures/.
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   MCCOOL             merged .mcool matrix (hic2cool output)
#   GENOME_BED12       gene annotation in BED12 format
#   HISTONE_BW_DIR     bigwig directory <mark>/<mark>_FE.bw
#   RNA_BW_CRO         RNA-seq bigwig, control (Cro)
#   RNA_BW_A485        RNA-seq bigwig, treated (10uM)
#   LOOPS_BEDPE        optional, default <HI_C_DIR>/loop/K1_merged.hiccups/merged_loops.bedpe
#   HIC_FIG_DIR        optional, default <HI_C_DIR>/figures
import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import trackc as tc

warnings.simplefilter("ignore")

HI_C_DIR = os.environ.get("HI_C_DIR", os.path.join(os.getcwd(), "hi_c"))
OUTDIR = os.environ.get("HIC_FIG_DIR", os.path.join(HI_C_DIR, "figures"))
os.makedirs(OUTDIR, exist_ok=True)

LOOPS = os.environ.get("LOOPS_BEDPE", os.path.join(HI_C_DIR, "loop/K1_merged.hiccups/merged_loops.bedpe"))
CRO_BW = os.environ.get("RNA_BW_CRO", "")
A485_BW = os.environ.get("RNA_BW_A485", "")
BWD = os.environ.get("HISTONE_BW_DIR", "")
GENE_BED12 = os.environ.get("GENOME_BED12", "")
K1_MCOOL = (os.environ.get("MCCOOL", "") + "::/resolutions/10000")
ONTAD_DIR = os.path.join(HI_C_DIR, "TAD/10kb.onTAD")

# (gene, chrom, start, end, rna_vmax, extra_genes_in_gene_track)
# extra_genes: genes to draw in the gene track in addition to the target gene
# (e.g. the distal partner of a P-P loop). P-P partners can also be added
# automatically via find_pp_distal_genes() when the loop-gene tables are present.
GENES = [
    ("LOC_00004064", "scaffold_11", 850_000, 1_755_000, 5000, ["LOC_00004056"]),
    ("LOC_00010466", "scaffold_3", 211_685_000, 212_925_000, 500, []),
    ("LOC_00014079", "scaffold_5", 85_780_000, 87_045_000, 800, []),
    ("LOC_00012676", "scaffold_4", 165_130_000, 166_210_000, 5000, ["LOC_00012660"]),
    # batch3 loci (README-batch3.md plot ranges)
    ("LOC_00010536", "scaffold_3", 229_471_055, 230_455_000, 5000, []),
    ("LOC_00014813", "scaffold_6", 1_150_000, 2_385_000, 1500, []),
    ("LOC_00002287", "scaffold_1", 301_230_000, 302_420_000, 3000, []),
    ("LOC_00001911", "scaffold_1", 262_838_036, 263_780_000, 4000, []),
    ("LOC_00017828", "scaffold_7", 111_520_000, 112_375_000, 5000, []),
]

# ---- P-P loop distal-gene detection (optional; needs loop-gene classification) ----
LOOP_GENE_DETAIL = os.path.join(HI_C_DIR, "loop/gene_classification/loop_gene_detail.tsv")
GENE_CLF = os.path.join(HI_C_DIR, "loop/gene_classification/gene_loop_classification.tsv")
G5K1 = os.path.join(HI_C_DIR, "loop/K1_merged.hiccups/anchor1.intersect.gene5kb.bed")
G5K2 = os.path.join(HI_C_DIR, "loop/K1_merged.hiccups/anchor2.intersect.gene5kb.bed")

_HAVE_PP = all(os.path.exists(f) for f in (G5K1, G5K2, LOOP_GENE_DETAIL, GENE_CLF))
if _HAVE_PP:
    _g5k1 = pd.read_csv(G5K1, sep="\t", header=None)
    _g5k1.columns = ["chr", "as", "ae", "pair", "gchr", "gs", "ge", "gene"]
    _g5k2 = pd.read_csv(G5K2, sep="\t", header=None)
    _g5k2.columns = ["chr", "as", "ae", "pair", "gchr", "gs", "ge", "gene"]
    G5K = {"anchor1": _g5k1, "anchor2": _g5k2}
    _lg = pd.read_csv(LOOP_GENE_DETAIL, sep="\t", dtype=str)
    _clf = pd.read_csv(GENE_CLF, sep="\t", dtype=str)
else:
    G5K = {}
    _lg = pd.DataFrame()
    _clf = pd.DataFrame()


def _anchor_hits_gene(side, pair, chr_, s, e):
    df = G5K[side]
    return len(df[(df["pair"] == pair) & (df["chr"] == chr_) & (df["as"] == s) & (df["ae"] == e)]) > 0


def find_pp_distal_genes(gene):
    """For loops of `gene` whose BOTH anchors hit a gene TSS+/-5kb (P-P type),
    return the gene IDs on the distal anchor."""
    out = []
    if not _HAVE_PP:
        return out
    cl = _clf[_clf["geneID"] == gene]
    if len(cl) == 0 or int(cl["loop_count"].iloc[0]) == 0:
        return out
    for pair in cl["loop_list"].iloc[0].split(","):
        for _, r in _lg[(_lg["loopPair"] == pair) & (_lg["geneID"] == gene)].iterrows():
            gs = r["anchorType"]
            other_side = r["distal_anchorType"]
            ga = (r["anchor_chr"], int(r["anchor_start"]), int(r["anchor_end"]))
            oth = (r["distal_anchor_chr"], int(r["distal_anchor_start"]), int(r["distal_anchor_end"]))
            if not (_anchor_hits_gene(gs, pair, *ga) and _anchor_hits_gene(other_side, pair, *oth)):
                continue
            df = G5K[other_side]
            sub = df[(df["pair"] == pair) & (df["chr"] == oth[0]) & (df["as"] == oth[1]) & (df["ae"] == oth[2])]
            for gid in sub["gene"].tolist():
                if gid != gene:
                    out.append(gid)
    return sorted(set(out))

LEVEL_COLORS = {1: "blue", 2: "#44C413", 3: "purple", 4: "#FE9900", 5: "#37ECF3"}
TC_MARKS = ["H3K27ac", "H3K4me3", "H3K4me1", "H3K36me3", "H3K27me3", "H3K9me3"]


# ---------------- state_trackc (loop links + histone marks + RNA + gene) ----------------
def make_state_trackc(gene, chrom, start, end, rna_vmax, extra_genes=()):
    extra_genes = list(extra_genes) + find_pp_distal_genes(gene)
    loops = pd.read_table(LOOPS, comment="#", header=None)
    gene_bed12 = pd.read_table(GENE_BED12, header=None)
    bws = {m: pyBigWig.open(f"{BWD}/{m}/{m}_FE.bw") for m in TC_MARKS}
    cro = pyBigWig.open(CRO_BW)
    a485 = pyBigWig.open(A485_BW)
    regions = [f"{chrom}:{start}-{end}"]

    ten = tc.tenon(figsize=(8, 0.6))
    ten.add(pos="bottom", height=1)
    for _ in range(8):
        ten.add(pos="bottom", height=1, hspace=0.1)
    ten.add(pos="bottom", height=0.8, hspace=0.1)
    ten.add(pos="bottom", height=0.5, hspace=0.2)

    tc.pl.links_track(ax=ten.axs(0), data=loops, label="loop", regions=regions,
                      color=["#66AC84", "tab:purple"], maxrange=3_000_000, anchor="inside")
    tc.pl.bw_track(bws["H3K27ac"], ten.axs(1), regions=regions, vmax=50,
                   label="H3K27ac", binsize=100, primary_col="tab:blue")
    tc.pl.bw_track(bws["H3K4me3"], ten.axs(2), regions=regions, vmax=50,
                   label="H3K4me3", binsize=100, primary_col="tab:blue")
    tc.pl.bw_track(bws["H3K4me1"], ten.axs(3), regions=regions, vmax=50,
                   label="H3K4me1", binsize=100, primary_col="tab:blue")
    tc.pl.bw_track(bws["H3K36me3"], ten.axs(4), regions=regions, vmax=50,
                   label="H3K36me3", binsize=100, primary_col="tab:blue")
    tc.pl.bw_track(bws["H3K27me3"], ten.axs(5), regions=regions, vmax=10,
                   label="H3K27me3", binsize=100, primary_col="tab:blue")
    tc.pl.bw_track(bws["H3K9me3"], ten.axs(6), regions=regions, vmax=10,
                   label="H3K9me3", binsize=100, primary_col="tab:blue")
    tc.pl.bw_track(a485, ten.axs(7), regions=regions, vmax=rna_vmax,
                   label="RNA-seq A485", binsize=100, primary_col="tab:orange")
    tc.pl.bw_track(cro, ten.axs(8), regions=regions, vmax=rna_vmax,
                   label="RNA-seq Cro", binsize=100, primary_col="tab:blue")
    _g = gene_bed12[gene_bed12[3].isin([gene] + list(extra_genes))].copy()
    tc.pl.gene_track(ax=ten.axs(9), bed12=_g, regions=regions, line=1,
                     ax_on=True, gene_fontszie=10)
    tc.pl.multi_scale_track(ten.axs(10), regions=regions, scale_adjust="Mb",
                            intervals=1, tick_fontsize=10)

    out = os.path.join(OUTDIR, f"Figure6_{gene}_state_trackc.pdf")
    tc.savefig(out)
    for b in list(bws.values()) + [cro, a485]:
        b.close()
    print("saved:", out)


# ---------------- markline (Hi-C triangle + onTAD domain levels) ----------------
def load_ontad_levels(chrom, lo, hi):
    bed = os.path.join(ONTAD_DIR, f"onTAD.K1.{chrom}.bed")
    levels = defaultdict(list)
    for l in open(bed):
        p = l.split()
        if not p or p[0].startswith("track") or p[0].startswith("#"):
            continue
        if p[0] != chrom:
            continue
        s, e, lv = int(p[1]), int(p[2]), p[-1]
        if e <= lo or s >= hi:
            continue
        s, e = max(s, lo), min(e, hi)
        levels[lv].append(f"{chrom}:{s}-{e}")
    return levels


def auto_maxrange(offdiag_pos):
    if len(offdiag_pos) == 0:
        return 100.0
    return float(np.percentile(offdiag_pos, 99))


def make_markline(gene, chrom, start, end):
    regions = [f"{chrom}:{start}-{end}"]
    levels = load_ontad_levels(chrom, start, end)
    print(f"{gene} onTAD levels:", {k: len(v) for k, v in levels.items()})

    z = tc.tl.extractContactRegions(clr=K1_MCOOL, row_regions=regions)
    mat = np.asarray(z.cmat, dtype=float)
    offdiag = mat[np.triu_indices_from(mat, k=1)]
    offdiag = offdiag[np.isfinite(offdiag)]
    offdiag = offdiag[offdiag > 0]
    maxrange = auto_maxrange(offdiag)
    print(f"{gene} auto maxrange = {maxrange:.1f}")

    ten = tc.tenon(figsize=(6, 1))
    ten.add(pos="bottom", height=3, hspace=0.05)
    tc.pl.mapC(ax=ten.axs(0), mat=z.cmat, map_type="triangle",
               maxrange=maxrange, minrange=0, ax_on=False, height=220)
    for lv in sorted(levels, key=lambda x: int(x[5:]) if x.startswith("level") else 999):
        n = int(lv[5:]) if lv.startswith("level") else lv
        for region in levels[lv]:
            tc.pl.mapc_markline(ax=ten.axs(0), row_regions=z.row_regions,
                                map_type="triangle", mark_regions=[region],
                                binsize=10000, map_order=0, symmetric=False,
                                show_regions_edge=False, linecolor=LEVEL_COLORS[n],
                                linewidth=1.5, linestyle="dotted")
    pdf = os.path.join(OUTDIR, f"Figure6_{gene}_markline.pdf")
    png = os.path.join(OUTDIR, f"Figure6_{gene}_markline.png")
    ten.fig.savefig(pdf, bbox_inches="tight")
    ten.fig.savefig(png, bbox_inches="tight", dpi=220)
    print("saved:", pdf)
    print("saved:", png)


if __name__ == "__main__":
    for gene, chrom, start, end, rna_vmax, extra_genes in GENES:
        print("=" * 40)
        print(f"gene={gene}  {chrom}:{start:,}-{end:,}  RNA vmax={rna_vmax}  extra_genes={extra_genes}")
        make_state_trackc(gene, chrom, start, end, rna_vmax, extra_genes)
        make_markline(gene, chrom, start, end)
    print("\nALL DONE ->", OUTDIR)