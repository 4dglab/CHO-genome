#!/usr/bin/env python3
# PURPOSE: mmseqs clustering + genomic overlap redundancy removal
# INPUT:   EviAnn/AGAT GFF
# OUTPUT:  deduplicated GFF
# TOOL:    mmseqs2


"""Redundancy removal: mmseqs clusters + genomic overlap rules (B4 step).

Rules (per approved plan):
  1. Same cluster + same strand + genomic overlap  -> keep stronger evidence, drop weaker
  2. Same cluster + opposite strand + overlap      -> keep both (antisense, never drop)
  3. Same-strand nested small gene (fully inside a larger gene span):
       - no independent transcript evidence (Evidence not in {complete, transcript_only})
         -> drop the nested gene
       - with independent transcript evidence -> keep (genuine nested gene)
  4. Same-strand CDS partial overlap (not nested), weaker evidence -> drop weaker
  lncRNA genes are always kept.

Evidence rank: complete(3) > transcript_only(2) > protein_only(1) > ab_initio/augustus(0).

Usage:
  overlap_analysis.py --gff <genomic.gff> --clusters <clusters.tsv> \
      --outdir <redundancy/>

Outputs (in outdir):
  dropped_redundant.tsv   gene_id  reason  evidence  partner  detail
  nested_kept.tsv         genuine nested genes kept (report only)
  antisense_kept.tsv      antisense overlapping pairs kept
  overlap_summary.tsv     counts by decision type
"""
import argparse
import os
import re
from collections import defaultdict


ATTR_RE = re.compile(r"([^=;]+)=([^;]+)")


def parse_attrs(attr_str):
    d = {}
    for k, v in ATTR_RE.findall(attr_str):
        d[k] = v
    return d


def load_gff(path):
    """Return (genes, transcripts) dicts.

    genes: gene_id -> dict(seqid,start,end,strand,source,biotype,evidence_rank,
                           evidence_tag,n_transcripts,max_cds_len,cds_list)
    transcripts: tx_id -> gene_id  (mRNA / transcript / lnc_RNA features)
    """
    genes = {}
    tx2gene = {}
    gene_order = []
    cds_by_tx = defaultdict(list)
    gene_features = {}  # keep raw gene attribute line for later (not used here)

    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        seqid, source, ftype, start, end, _, strand, _, attr = f
        start, end = int(start), int(end)
        a = parse_attrs(attr)
        if ftype == "gene":
            gid = a.get("ID")
            if gid is None:
                continue
            biotype = a.get("gene_biotype", "protein_coding")
            genes[gid] = {
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "source": source,
                "biotype": biotype,
                "evidence_rank": 0,
                "evidence_tag": "none",
                "n_transcripts": 0,
                "max_cds_len": 0,
            }
            gene_order.append(gid)
        elif ftype in ("mRNA", "transcript", "lnc_RNA"):
            tid = a.get("ID")
            parent = a.get("Parent")
            if tid is None or parent is None:
                continue
            tx2gene[tid] = parent
            if parent in genes:
                genes[parent]["n_transcripts"] += 1
                ev = a.get("Evidence", "none")
                rank = {"complete": 3, "transcript_only": 2, "protein_only": 1}.get(ev, 0)
                if rank > genes[parent]["evidence_rank"]:
                    genes[parent]["evidence_rank"] = rank
                    genes[parent]["evidence_tag"] = ev
        elif ftype == "CDS":
            tid = a.get("Parent")
            if tid:
                cds_by_tx[tid].append((start, end))
    for tid, spans in cds_by_tx.items():
        if tid in tx2gene:
            gid = tx2gene[tid]
            clen = sum(e - s + 1 for s, e in spans)
            if clen > genes[gid]["max_cds_len"]:
                genes[gid]["max_cds_len"] = clen
    return genes, tx2gene, gene_order


def load_clusters(path):
    """Return dict gene_id -> cluster_id using protein(tx) names -> gene mapping."""
    clusters = {}
    return clusters


def load_cluster_groups(path, tx2gene, genes):
    """Map gene_id -> cluster_id. Skip members whose tx not in tx2gene or gene is lncRNA."""
    gene2cluster = {}
    cluster_members = defaultdict(list)
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        rep, member = parts[0], parts[1]
        g = tx2gene.get(member)
        if g is None or g not in genes:
            continue
        if genes[g]["biotype"] == "lncRNA":
            continue
        gene2cluster[g] = rep
        cluster_members[rep].append(g)
    return gene2cluster, cluster_members


def overlap(a, b):
    return a["start"] <= b["end"] and b["start"] <= a["end"]


def cds_overlap_bp(a, b, cds_by_gene):
    spans_a = cds_by_gene.get(a, [])
    spans_b = cds_by_gene.get(b, [])
    tot = 0
    for s1, e1 in spans_a:
        for s2, e2 in spans_b:
            ov = min(e1, e2) - max(s1, s2) + 1
            if ov > 0:
                tot += ov
    return tot


def is_weaker(a, b, genes):
    ra = genes[a]["evidence_rank"]
    rb = genes[b]["evidence_rank"]
    if ra != rb:
        return ra < rb
    return genes[a]["max_cds_len"] < genes[b]["max_cds_len"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    genes, tx2gene, gene_order = load_gff(args.gff)
    gene2cluster, cluster_members = load_cluster_groups(args.clusters, tx2gene, genes)

    # CDS spans per gene (union of all transcripts)
    cds_by_gene = defaultdict(list)
    for line in open(args.gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = parse_attrs(f[8])
        tid = a.get("Parent")
        if tid in tx2gene:
            cds_by_gene[tx2gene[tid]].append((int(f[3]), int(f[4])))

    dropped = []      # (gene, reason, partner, detail)
    kept = set(gene_order)
    nested_kept = []
    antisense_kept = []

    # ---- Pass 1: mmseqs cluster redundancy (same-strand overlapping) ----
    for rep, members in cluster_members.items():
        if len(members) < 2:
            continue
        for i in range(len(members)):
            gi = members[i]
            if gi not in kept:
                continue
            for j in range(i + 1, len(members)):
                gj = members[j]
                if gj not in kept:
                    continue
                a, b = genes[gi], genes[gj]
                if a["seqid"] != b["seqid"] or not overlap(a, b):
                    continue
                if a["strand"] == b["strand"]:
                    weaker = gi if is_weaker(gi, gj, genes) else gj
                    stronger = gj if weaker == gi else gi
                    dropped.append((weaker, "cluster_redundant",
                                    stronger,
                                    "mmseqs cluster %s, %s-%s bp, overlap" % (
                                        rep, genes[weaker]["max_cds_len"],
                                        genes[stronger]["max_cds_len"])))
                    kept.discard(weaker)
                else:
                    antisense_kept.append((gi, gj))

    # ---- Pass 2: nested genes (same strand, fully contained) ----
    same_strand = defaultdict(list)
    for g in gene_order:
        if g in kept:
            same_strand[(genes[g]["seqid"], genes[g]["strand"])].append(g)
    for key, gs in same_strand.items():
        gs_sorted = sorted(gs, key=lambda g: (genes[g]["start"], -(genes[g]["end"])))
        for i, outer in enumerate(gs_sorted):
            if outer not in kept:
                continue
            o = genes[outer]
            for inner in gs_sorted[i + 1:]:
                if inner not in kept:
                    continue
                n = genes[inner]
                # fully inside outer span
                if n["start"] >= o["start"] and n["end"] <= o["end"]:
                    # if lncRNA inner: keep
                    if genes[inner]["biotype"] == "lncRNA":
                        continue
                    indep = genes[inner]["evidence_tag"] in ("complete", "transcript_only")
                    if not indep:
                        dropped.append((inner, "nested_no_evidence", outer,
                                        "inside %s (%s-%s)" % (
                                            outer, o["start"], o["end"])))
                        kept.discard(inner)
                    else:
                        nested_kept.append((inner, outer))

    # ---- Pass 3: same-strand CDS partial overlap (weaker evidence) ----
    for key, gs in same_strand.items():
        gs_sorted = sorted(gs, key=lambda g: genes[g]["start"])
        for i in range(len(gs_sorted)):
            gi = gs_sorted[i]
            if gi not in kept:
                continue
            for j in range(i + 1, len(gs_sorted)):
                gj = gs_sorted[j]
                if gj not in kept:
                    continue
                a, b = genes[gi], genes[gj]
                if a["end"] < b["start"]:
                    break
                if a["start"] <= b["start"] and a["end"] >= b["end"]:
                    continue  # handled by nested pass
                ov = cds_overlap_bp(gi, gj, cds_by_gene)
                if ov == 0:
                    continue
                min_cds = min(genes[gi]["max_cds_len"], genes[gj]["max_cds_len"])
                if min_cds > 0 and ov / min_cds >= 0.5 and \
                        genes[gi]["evidence_rank"] != genes[gj]["evidence_rank"]:
                    weaker = gi if is_weaker(gi, gj, genes) else gj
                    stronger = gj if weaker == gi else gi
                    dropped.append((weaker, "cds_partial_overlap", stronger,
                                    "CDS overlap %d bp (>=50%% of shorter)" % ov))
                    kept.discard(weaker)

    # ---- Output ----
    with open(os.path.join(args.outdir, "dropped_redundant.tsv"), "w") as out:
        out.write("gene_id\tseqid\tstart\tend\tstrand\tsource\tbiotype\tevidence\treason\tpartner\tdetail\n")
        for g, reason, partner, detail in dropped:
            ge = genes[g]
            out.write("\t".join(map(str, [
                g, ge["seqid"], ge["start"], ge["end"], ge["strand"],
                ge["source"], ge["biotype"], ge["evidence_tag"], reason, partner, detail])) + "\n")

    with open(os.path.join(args.outdir, "nested_kept.tsv"), "w") as out:
        out.write("inner\touter\n")
        for inner, outer in nested_kept:
            out.write("%s\t%s\n" % (inner, outer))

    with open(os.path.join(args.outdir, "antisense_kept.tsv"), "w") as out:
        out.write("gene_a\tgene_b\n")
        for a, b in antisense_kept:
            out.write("%s\t%s\n" % (a, b))

    with open(os.path.join(args.outdir, "kept_genes.tsv"), "w") as out:
        for g in gene_order:
            out.write(g + "\tkept" if g in kept else g + "\tdropped")
            out.write("\n")

    summary = {
        "total_genes": len(gene_order),
        "cluster_redundant_dropped": sum(1 for d in dropped if d[1] == "cluster_redundant"),
        "nested_no_evidence_dropped": sum(1 for d in dropped if d[1] == "nested_no_evidence"),
        "cds_partial_overlap_dropped": sum(1 for d in dropped if d[1] == "cds_partial_overlap"),
        "total_dropped": len(dropped),
        "total_kept": len(kept),
        "genuine_nested_kept": len(nested_kept),
        "antisense_overlap_kept": len(antisense_kept),
        "genes_in_clusters": len(gene2cluster),
    }
    with open(os.path.join(args.outdir, "overlap_summary.tsv"), "w") as out:
        for k, v in summary.items():
            out.write("%s\t%s\n" % (k, v))

    print("OVERLAP_ANALYSIS_DONE")
    for k, v in summary.items():
        print("  %-30s %s" % (k, v))


if __name__ == "__main__":
    main()
