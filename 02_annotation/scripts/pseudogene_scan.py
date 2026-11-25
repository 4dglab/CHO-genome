#!/usr/bin/env python3
# PURPOSE: Pseudogene / truncated model filtering
# INPUT:   deduplicated GFF
# OUTPUT:  filtered GFF
# TOOL:    custom


"""Pseudogene / truncated-model filtering (B4 step).

Signals:
  - internal stop codon (in-frame '*')            -> pseudogene   -> DROP
  - broken ORF (CDS length not divisible by 3)     -> pseudogene   -> DROP
  - missing start AND stop (translated ORF)        -> truncated    -> DROP
  - EviAnn pseudo_label (gene_biotype=pseudogene)  -> pseudogene   -> DROP

Completeness logic:
  - EviAnn transcripts carry explicit StartCodon/StopCodon attributes in the GFF
    -> treated as complete unless an internal stop is found.
  - Transcripts without those attributes (e.g. Augustus ab initio) are translated
    from the genome; a proper ORF must start with M and end with a stop codon
    (TAA/TAG/TGA). Missing start only  -> partial (keep, flagged);
    missing stop only -> incomplete (keep, flagged); missing both -> truncated (drop).

Usage:
  pseudogene_scan.py --proteins <primary.proteins.fa> --gff <genomic.gff> \
      --genome <softmasked.fa> --pseudo-label <pseudo_label.gff> \
      --outdir <pseudogene/>

Outputs (in outdir):
  pseudogene_dropped.tsv    gene_id reason detail
  partial_flagged.tsv       kept but flagged partial/incomplete genes
  pseudogene_summary.tsv    counts
"""
import argparse
import os
import re
from collections import defaultdict

GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
COMP = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
STOPS = ("TAA", "TAG", "TGA")


def load_genome(path):
    seqs = {}
    cur = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def revcomp(s):
    return "".join(COMP[b] for b in reversed(s))


def translate(cds):
    aa = []
    for i in range(0, len(cds) - 2, 3):
        aa.append(GENETIC_CODE.get(cds[i:i+3], "X"))
    return "".join(aa)


def parse_fasta(path):
    seqs = {}
    cur = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--pseudo-label", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    genome = load_genome(args.genome)
    proteins = parse_fasta(args.proteins)

    # transcripts -> gene, seqid, strand, CDS exon list (ordered), GFF start/stop attrs
    tx_info = {}
    cds_by_tx = defaultdict(list)
    for line in open(args.gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = dict(re.findall(r"([^=;]+)=([^;]+)", f[8]))
        if f[2] in ("mRNA", "transcript", "lnc_RNA"):
            if a.get("ID") and a.get("Parent"):
                tx_info[a["ID"]] = {
                    "gene": a["Parent"], "seqid": f[0], "strand": f[6],
                    "start_codon": a.get("StartCodon", ""),
                    "stop_codon": a.get("StopCodon", ""),
                }
        elif f[2] == "CDS" and a.get("Parent") in tx_info:
            cds_by_tx[a["Parent"]].append((int(f[3]), int(f[4])))

    # EviAnn pseudo_label genes
    pseudo_label = set()
    for line in open(args.pseudo_label):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(re.findall(r"([^=;]+)=([^;]+)", f[8]))
        if "pseudogene" in a.get("gene_biotype", "") and a.get("ID"):
            pseudo_label.add(a["ID"])

    dropped = []
    flagged = []
    gff_attrs_complete = 0
    translated_ok = 0
    translated_nostart_nostop = 0

    for tx, info in tx_info.items():
        gid = info["gene"]
        exons = sorted(cds_by_tx.get(tx, []), key=lambda x: (x[0], x[1]))
        if not exons:
            continue
        # concatenate CDS (5'->3' on strand)
        seq = genome[info["seqid"]]
        cds = ""
        if info["strand"] == "+":
            cds = "".join(seq[s - 1:e] for s, e in exons)
        else:
            cds = "".join(revcomp(seq[s - 1:e]) for s, e in reversed(exons))
        if len(cds) % 3 != 0:
            dropped.append((gid, "broken_orf",
                            "CDS length %d not divisible by 3 (%s)" % (len(cds), tx)))
            continue
        aa = translate(cds)
        # translate() includes the terminal stop codon as '*'
        has_stop = aa.endswith("*")
        body = aa[:-1] if has_stop else aa
        if "*" in body:
            dropped.append((gid, "internal_stop",
                            "in-frame stop codon in %s (CDS %d bp)" % (tx, len(cds))))
            continue
        # completeness: EviAnn explicitly annotated start/stop codons
        if info["start_codon"] and info["stop_codon"]:
            gff_attrs_complete += 1
            continue  # complete by annotation
        # otherwise translate-based check (Augustus etc.)
        translated_ok += 1
        has_start = body.startswith("M")
        if not has_start and not has_stop:
            translated_nostart_nostop += 1
            dropped.append((gid, "truncated_no_start_stop",
                            "%s: no start (M) and no stop codon" % tx))
        elif not has_start:
            flagged.append((gid, "no_start_codon",
                            "%s: no start methionine (has stop)" % tx))
        elif not has_stop:
            flagged.append((gid, "no_stop_codon",
                            "%s: has start but no stop codon" % tx))

    for gid in pseudo_label:
        dropped.append((gid, "eviann_pseudo_label", "marked as processed pseudogene"))

    with open(os.path.join(args.outdir, "pseudogene_dropped.tsv"), "w") as out:
        out.write("gene_id\treason\tdetail\n")
        for g, r, d in dropped:
            out.write("%s\t%s\t%s\n" % (g, r, d))

    with open(os.path.join(args.outdir, "partial_flagged.tsv"), "w") as out:
        out.write("gene_id\treason\tdetail\n")
        for g, r, d in flagged:
            out.write("%s\t%s\t%s\n" % (g, r, d))

    summary = {
        "transcripts_analyzed": len(tx_info),
        "complete_by_attrs": gff_attrs_complete,
        "translation_checked": translated_ok,
        "internal_stop_pseudogene": sum(1 for d in dropped if d[1] == "internal_stop"),
        "broken_orf": sum(1 for d in dropped if d[1] == "broken_orf"),
        "truncated_no_start_stop": sum(1 for d in dropped if d[1] == "truncated_no_start_stop"),
        "eviann_pseudo_label": sum(1 for d in dropped if d[1] == "eviann_pseudo_label"),
        "pseudogene_dropped": len(dropped),
        "partial_no_start_kept": sum(1 for d in flagged if d[1] == "no_start_codon"),
        "incomplete_no_stop_kept": sum(1 for d in flagged if d[1] == "no_stop_codon"),
    }
    with open(os.path.join(args.outdir, "pseudogene_summary.tsv"), "w") as out:
        for k, v in summary.items():
            out.write("%s\t%s\n" % (k, v))

    print("PSEUDOGENE_SCAN_DONE")
    for k, v in summary.items():
        print("  %-30s %s" % (k, v))


if __name__ == "__main__":
    main()
