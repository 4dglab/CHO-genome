#!/usr/bin/env python3
# PURPOSE: PUR-associated gene classification
# INPUT:   PUR BED + gene GFF
# OUTPUT:  pur_gene_classification.tsv
# TOOL:    custom

"""Classify protein-coding genes located in PUR (previously unassembled regions).

PUR membership: a gene is considered "in PUR" if >= 75% of its genomic span
overlaps PUR intervals (pur_full.ge1kb.bed, 0-based half-open).

Gene classification (final NR-threshold version):
  previously annotated      : protein-coding gene NOT in the 785 novel set
                              (= 20139 - 785 = 19354, matched to reference)
  unannotated               : novel protein-coding gene with a NR hit passing
                              thresholds (category "unannotated")
  previously unidentified   : novel protein-coding gene with no NR hit or a
                              hit below threshold (categories "previously_
                              unidentified" + "nonNRhit")

Inputs:
  --pur-bed     PUR intervals (BED, 0-based half-open)
  --gff         high-confidence GFF (protein-coding incl L4, 1-based closed)
  --class       novel_protein_coding_nr_classification.tsv (785 genes)
Outputs:
  --out-summary summary TSV (with explanation header)
  --out-detail  per-gene TSV
"""
import argparse
import bisect
import collections

ap = argparse.ArgumentParser()
ap.add_argument("--pur-bed", required=True)
ap.add_argument("--gff", required=True)
ap.add_argument("--class", dest="cls", required=True)
ap.add_argument("--out-summary", required=True)
ap.add_argument("--out-detail", required=True)
ap.add_argument("--frac", type=float, default=0.75, help="min gene-span overlap fraction (default 0.75)")
a = ap.parse_args()

PUR_FRAC = a.frac

# ---------- PUR intervals (0-based half-open), sorted per scaffold ----------
pur = collections.defaultdict(list)
for line in open(a.pur_bed):
    p = line.rstrip("\n").split("\t")
    if len(p) < 3:
        continue
    pur[p[0]].append((int(p[1]), int(p[2])))
pur_starts = {k: sorted(x[0] for x in v) for k, v in pur.items()}
for k in pur:
    pur[k].sort()


def pur_overlap(seqid, gs, ge):
    """Overlap (bp) of 0-based half-open interval [gs, ge) with PUR."""
    iv = pur.get(seqid)
    if not iv:
        return 0
    i = bisect.bisect_right(pur_starts[seqid], gs)
    j = bisect.bisect_left(pur_starts[seqid], ge)
    total = 0
    for s, e in iv[max(0, i - 1):j]:
        total += max(0, min(e, ge) - max(s, gs))
    return total


# ---------- gene classification ----------
novel_cat = {}
for line in open(a.cls):
    if line.startswith("gene_id"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) >= 2:
        novel_cat[f[0]] = f[1]


def classify(gid):
    c = novel_cat.get(gid)
    if c is None:
        return "previously_annotated"
    if c == "unannotated":
        return "unannotated"
    return "previously_unidentified"  # previously_unidentified + nonNRhit


# ---------- iterate genes ----------
genes = []
n_pc = 0
for line in open(a.gff):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] != "gene":
        continue
    at = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
    bio = at.get("gene_biotype", "L4")
    if not (bio == "protein_coding" or bio == "L4"):
        continue
    gid = at.get("ID")
    if not gid:
        continue
    n_pc += 1
    s, e = int(p[3]), int(p[4])
    length = e - s + 1
    ov = pur_overlap(p[0], s - 1, e)  # convert 1-based closed -> 0-based half-open
    frac = ov / length
    genes.append({
        "gene_id": gid, "seqid": p[0], "start": s, "end": e,
        "length": length, "category": classify(gid),
        "pur_overlap": ov, "pur_frac": frac, "in_pur": frac >= PUR_FRAC,
    })

# ---------- summary ----------
cnt_all = collections.Counter(g["category"] for g in genes)
in_pur = [g for g in genes if g["in_pur"]]
cnt_pur = collections.Counter(g["category"] for g in in_pur)
total_pur = len(in_pur)

lines = []
lines.append("# protein-coding gene classification within PURs")
lines.append("# =====================================================================")
lines.append("# PUR definition: pur_full.ge1kb.bed (0-based half-open), genome-wide >=1 kb")
lines.append("# gene in PUR: gene span >= %.0f%% overlaps a PUR interval" % (PUR_FRAC * 100))
lines.append("# categories (final NR-threshold version):")
lines.append("#   previously annotated   : not among the 785 novel protein-coding genes")
lines.append("#                             protein-coding genes (matched to reference, %d - %d = %d)" % (
    n_pc, sum(cnt_all.get(x, 0) for x in ("unannotated", "previously_unidentified")),
    cnt_all["previously_annotated"]))
lines.append("#   unannotated            : novel protein-coding genes with an NR hit above threshold")
lines.append("#                             (the 'unannotated' class in")
lines.append("#                              novel_protein_coding_nr_classification.tsv)")
lines.append("#   previously unidentified: novel protein-coding genes without an NR hit")
lines.append("#                             (previously_unidentified + nonNRhit classes)")
lines.append("# denominator = total protein-coding genes within PURs")
lines.append("# =====================================================================")
lines.append("category\tcount_in_pur\tpct_of_pur\tcount_all")
for c in ("previously_annotated", "unannotated", "previously_unidentified"):
    lines.append("%s\t%d\t%.2f\t%d" % (c, cnt_pur[c], 100.0 * cnt_pur[c] / total_pur, cnt_all[c]))
lines.append("TOTAL\t%d\t100.00\t%d" % (total_pur, n_pc))
with open(a.out_summary, "w") as out:
    out.write("\n".join(lines) + "\n")

with open(a.out_detail, "w") as out:
    out.write("gene_id\tseqid\tstart\tend\tlength\tcategory\tpur_overlap_bp\tpur_frac\tin_pur_ge75\n")
    for g in sorted(genes, key=lambda x: (x["seqid"], x["start"])):
        out.write("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%.4f\t%s\n" % (
            g["gene_id"], g["seqid"], g["start"], g["end"], g["length"],
            g["category"], g["pur_overlap"], g["pur_frac"],
            "YES" if g["in_pur"] else "no"))

print("protein-coding genes total:", n_pc)
print("all-category counts:", dict(cnt_all))
print("in-PUR total:", total_pur)
print("in-PUR by category:", dict(cnt_pur))
for c in ("previously_annotated", "unannotated", "previously_unidentified"):
    print("  %-24s %4d  %.2f%%" % (c, cnt_pur[c], 100.0 * cnt_pur[c] / total_pur))
print("wrote:", a.out_summary)
print("wrote:", a.out_detail)
