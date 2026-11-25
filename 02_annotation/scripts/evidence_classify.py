#!/usr/bin/env python3
# PURPOSE: Evidence tier classification L1-L4 (Iso-seq/RNA-seq/homology)
# INPUT:   filtered GFF + BAM
# OUTPUT:  evidence-tiered GFF
# TOOL:    custom


"""Transcript evidence classification L1-L4 (B4 step).

Evidence levels (for the manuscript table):
  L1  Iso-seq + RNA-seq + protein            (EviAnn Evidence=complete AND iso-seq overlap)
  L2  Transcriptome support (RNA-seq/Iso-seq) (Evidence=complete w/o iso-seq, or transcript_only)
  L3  Homology protein only                  (Evidence=protein_only)
  L4  ab initio only                         (Augustus source)
  lncRNA                                      (separate category, kept)

Iso-seq support = any Iso-seq read (pre-aligned BAM) overlapping the gene span,
detected via bedtools intersect (read count per gene).

Usage:
  evidence_classify.py --gff <genomic.gff> --isoseq-bam <isoseq_aligned.bam> \
      --bedtools <path-to-bedtools> --outdir <evidence/>

Outputs (in outdir):
  evidence_levels.tsv       gene-level table
  evidence_summary.tsv      counts per level
  genomic.evidence.gff      full GFF with Evidence_level appended to gene lines
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict, Counter


ATTR_RE = re.compile(r"([^=;]+)=([^;]+)")


def parse_attrs(a):
    return dict(ATTR_RE.findall(a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--isoseq-bam", required=True)
    ap.add_argument("--bedtools", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if not args.bedtools:
        args.bedtools = os.environ.get("CONDA_BASE", "/path/to/anaconda3") + "/envs/HiTE/bin/bedtools"
    if not os.path.exists(args.bedtools):
        sys.exit("bedtools not found: %s" % args.bedtools)

    # ---- collect gene spans + evidence from GFF ----
    genes = {}
    gene_order = []
    tx2gene = {}
    exon_counts = Counter()
    for line in open(args.gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        seqid, source, ftype, start, end = f[0], f[1], f[2], int(f[3]), int(f[4])
        a = parse_attrs(f[8])
        if ftype == "gene":
            gid = a.get("ID")
            if gid:
                genes[gid] = {
                    "seqid": seqid, "start": start, "end": end,
                    "strand": f[6], "source": source,
                    "biotype": a.get("gene_biotype", "protein_coding"),
                    "evidence": "none", "evidence_rank": 0,
                }
                gene_order.append(gid)
        elif ftype in ("mRNA", "transcript", "lnc_RNA"):
            tid, parent = a.get("ID"), a.get("Parent")
            if tid and parent:
                tx2gene[tid] = parent
                ev = a.get("Evidence", "none")
                rank = {"complete": 3, "transcript_only": 2, "protein_only": 1}.get(ev, 0)
                if parent in genes and rank > genes[parent]["evidence_rank"]:
                    genes[parent]["evidence_rank"] = rank
                    genes[parent]["evidence"] = ev
        elif ftype == "exon":
            if a.get("Parent") in tx2gene:
                exon_counts[tx2gene[a["Parent"]]] += 1

    # ---- iso-seq overlap via bedtools ----
    bed_lines = []
    for gid, g in genes.items():
        bed_lines.append("%s\t%d\t%d\t%s" % (g["seqid"], g["start"] - 1, g["end"], gid))
    isoseq_counts = {}
    with tempfile.TemporaryDirectory() as tmp:
        bed = os.path.join(tmp, "genes.bed")
        with open(bed, "w") as out:
            out.write("\n".join(bed_lines) + "\n")
        cmd = [args.bedtools, "intersect", "-a", bed, "-b", args.isoseq_bam,
               "-wa", "-c"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("bedtools intersect failed:\n" + r.stderr)
        for line in r.stdout.splitlines():
            p = line.split("\t")
            if len(p) == 5:
                isoseq_counts[p[3]] = int(p[4])

    # ---- assign levels ----
    def assign(g):
        if g["source"] == "augustus":
            return "L4"
        if g["biotype"] == "lncRNA":
            return "lncRNA"
        ev = g["evidence"]
        iso = isoseq_counts.get(gid, 0) > 0
        if ev == "complete":
            return "L1" if iso else "L2"
        if ev == "transcript_only":
            return "L2"
        if ev == "protein_only":
            return "L3"
        return "L4"

    rows = []
    for gid in gene_order:
        g = genes[gid]
        level = assign(g)
        rows.append({
            "gene_id": gid, "seqid": g["seqid"], "start": g["start"], "end": g["end"],
            "strand": g["strand"], "source": g["source"], "biotype": g["biotype"],
            "evidence": g["evidence"], "isoseq_reads": isoseq_counts.get(gid, 0),
            "exons": exon_counts.get(gid, 0), "level": level,
        })

    # ---- outputs ----
    with open(os.path.join(args.outdir, "evidence_levels.tsv"), "w") as out:
        out.write("gene_id\tseqid\tstart\tend\tstrand\tsource\tbiotype\tevidence\tisoseq_reads\texons\tlevel\n")
        for r in rows:
            out.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n" % (
                r["gene_id"], r["seqid"], r["start"], r["end"], r["strand"],
                r["source"], r["biotype"], r["evidence"], r["isoseq_reads"],
                r["exons"], r["level"]))

    lvl2gid = defaultdict(list)
    gene2level = {}
    for r in rows:
        lvl2gid[r["level"]].append(r["gene_id"])
        gene2level[r["gene_id"]] = r["level"]
    with open(os.path.join(args.outdir, "evidence_summary.tsv"), "w") as out:
        out.write("level\tcount\n")
        for lvl in ("L1", "L2", "L3", "L4", "lncRNA"):
            out.write("%s\t%d\n" % (lvl, len(lvl2gid.get(lvl, []))))

    # ---- annotated GFF (add Evidence_level to gene lines) ----
    out_gff = os.path.join(args.outdir, "genomic.evidence.gff")
    with open(args.gff) as fin, open(out_gff, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                fout.write(line)
                continue
            a = parse_attrs(f[8])
            gid = a.get("ID")
            if gid and gid in gene2level:
                lvl = gene2level[gid]
                if "Evidence_level" not in f[8]:
                    f[8] = "%s;Evidence_level=%s" % (f[8], lvl)
            fout.write("\t".join(f) + "\n")

    print("EVIDENCE_CLASSIFY_DONE")
    for lvl in ("L1", "L2", "L3", "L4", "lncRNA"):
        print("  %-8s %d genes" % (lvl, len(lvl2gid.get(lvl, []))))


if __name__ == "__main__":
    main()
