#!/usr/bin/env python3
# PURPOSE: StringTie+TransDecoder CDS remap + Augustus training set
# INPUT:   merged transcripts
# OUTPUT:  training set
# TOOL:    custom


"""Build an Augustus training/test set.

Maps TransDecoder CDS (transcript coordinates, seqid = transcript name) back
onto genomic coordinates using the StringTie merged GFF3 exon structure, then
applies quality filters.

Filters:
- keep only complete ORFs (Name ~ "type:complete")
- CDS total length >= min_cds_len (300)
- drop genes whose CDS overlaps repeat regions >= repeat_frac (0.5)
- remove redundant gene pairs (reciprocal CDS overlap >= 0.8, keep longer)

Output: clean gene/mRNA/CDS GFF3 in GENOMIC coordinates (renumbered IDs,
##gff-version 3), split into train/test with a fixed seed.

Usage:
    augustus_build_training.py stringtie_merged.gff3 transdecoder.gff3 repeats.gff3 train.gff3 test.gff3
"""

import random
import re
import sys
from collections import defaultdict

MIN_CDS = 300
REPEAT_FRAC = 0.5
REDUNDANCY = 0.8


def parse_gff3(path):
    feats = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)
            feats.append(
                {
                    "seqid": cols[0],
                    "type": cols[2],
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "strand": cols[6],
                    "phase": cols[7],
                    "attrs": attrs,
                }
            )
    return feats


def parse_repeat_intervals(path):
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            intervals[cols[0]].append((int(cols[3]), int(cols[4])))
    for k in intervals:
        ivs = sorted(intervals[k])
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        intervals[k] = merged
    return intervals


def rep_frac(cds, rep_ivs):
    if not rep_ivs:
        return 0.0
    total = 0
    for s, e in cds:
        for rs, re in rep_ivs:
            total += max(0, min(e, re) - max(s, rs) + 1)
    cds_len = sum(e - s + 1 for s, e in cds)
    return total / cds_len if cds_len else 1.0


def map_transcript_to_genome(t_coords, exons_5to3, strand, seqid):
    """Map a transcript-coordinate interval [t1, t2] (1-based, 5'->3') to genomic coords.

    exons_5to3: list of (g_start, g_end) in 5'->3' order (minus strand: rightmost first).
    Returns genomic (start, end) with start <= end (strand carried separately).
    """
    t1, t2 = t_coords
    if t1 > t2:
        t1, t2 = t2, t1
    gpos = []
    cum = 0
    for gs, ge in exons_5to3:
        length = ge - gs + 1
        for t in range(max(t1, cum + 1), min(t2, cum + length) + 1):
            if strand == "+":
                gpos.append(gs + (t - cum - 1))
            else:
                gpos.append(ge - (t - cum - 1))
        cum += length
        if cum >= t2:
            break
    return min(gpos), max(gpos)


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    merged_gff, td_gff, rep_gff, out_train, out_test = sys.argv[1:]

    merged_feats = parse_gff3(merged_gff)
    strand_of = {}
    seqid_of = {}
    exons_5to3 = defaultdict(list)  # transcript_id -> [(g_start,g_end), ...] in 5'->3'
    for f in merged_feats:
        if f["type"] == "mRNA":
            strand_of[f["attrs"]["ID"]] = f["strand"]
            seqid_of[f["attrs"]["ID"]] = f["seqid"]
        elif f["type"] == "exon":
            pid = f["attrs"]["Parent"]
            if pid not in strand_of:
                strand_of[pid] = f["strand"]
            if pid not in seqid_of:
                seqid_of[pid] = f["seqid"]
            exons_5to3[pid].append((f["start"], f["end"]))
    for pid in exons_5to3:
        strand = strand_of.get(pid, "+")
        exons_5to3[pid].sort(key=lambda x: x[0], reverse=(strand == "-"))

    td_feats = parse_gff3(td_gff)
    genes = {}
    mrna_parent = {}
    cds_children = defaultdict(list)
    for f in td_feats:
        if f["type"] == "gene":
            genes[f["attrs"]["ID"]] = f
        elif f["type"] == "mRNA":
            mrna_parent[f["attrs"]["ID"]] = f["attrs"]["Parent"]
        elif f["type"] == "CDS":
            cds_children[f["attrs"]["Parent"]].append(f)

    rep_ivs = parse_repeat_intervals(rep_gff)

    selected = []  # (seqid, strand, gene_name, [(g_start,g_end),...])
    for gid, gf in genes.items():
        name = gf["attrs"].get("Name", "")
        if "type:complete" not in name:
            continue
        mrnas = [m for m, p in mrna_parent.items() if p == gid]
        if not mrnas:
            continue
        tid = mrnas[0]
        tname = re.sub(r"\.p\d+$", "", tid)
        if tname not in exons_5to3:
            continue
        cds = cds_children[tid]
        if not cds:
            continue
        strand = strand_of.get(tname, "+")
        if strand not in ("+", "-"):
            continue
        seqid = seqid_of.get(tname, "")
        if not seqid:
            continue
        g_cds = []
        for c in sorted(cds, key=lambda x: x["start"]):
            gs, ge = map_transcript_to_genome(
                (c["start"], c["end"]), exons_5to3[tname], strand, seqid
            )
            g_cds.append((gs, ge))
        cds_len = sum(e - s + 1 for s, e in g_cds)
        if cds_len < MIN_CDS:
            continue
        if rep_frac(g_cds, rep_ivs.get(seqid, [])) >= REPEAT_FRAC:
            continue
        mrna_feat = next((f for f in td_feats if f["type"] == "mRNA" and f["attrs"]["ID"] == tid), None)
        selected.append((seqid, strand, tid, g_cds, mrna_feat))

    selected.sort(key=lambda x: (x[0], x[3][0][0]))

    kept = []
    for seqid, strand, tid, cds, mrna_feat in selected:
        cds_len = sum(e - s + 1 for s, e in cds)
        dup = False
        removals = []
        for i, (seqid2, strand2, tid2, cds2, _) in enumerate(kept):
            if seqid != seqid2 or strand != strand2:
                continue
            cds2_len = sum(e - s + 1 for s, e in cds2)
            inter = sum(
                max(0, min(e, e2) - max(s, s2) + 1)
                for s, e in cds
                for s2, e2 in cds2
            )
            if inter / cds_len >= REDUNDANCY and inter / cds2_len >= REDUNDANCY:
                if cds_len <= cds2_len:
                    dup = True
                    break
                removals.append(i)
        if dup:
            continue
        kept = [k for i, k in enumerate(kept) if i not in removals]
        kept.append((seqid, strand, tid, cds, mrna_feat))

    random.seed(42)
    items = list(kept)
    random.shuffle(items)
    n_test = max(1, int(len(items) * 0.1))
    test_set = items[:n_test]
    train_set = items[n_test:]
    test_set.sort(key=lambda x: (x[0], x[3][0][0]))
    train_set.sort(key=lambda x: (x[0], x[3][0][0]))

    def write_gff3(path, genes, start_idx):
        with open(path, "w") as out:
            out.write("##gff-version 3\n")
            for idx, (seqid, strand, tid, cds, mrna_feat) in enumerate(genes, start_idx):
                gid = f"g{idx:06d}"
                g_start = min(s for s, e in cds)
                g_end = max(e for s, e in cds)
                out.write(f"{seqid}\tAugustusTrain\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\tID={gid}\n")
                out.write(f"{seqid}\tAugustusTrain\tmRNA\t{g_start}\t{g_end}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}\n")
                for s, e in sorted(cds):
                    out.write(f"{seqid}\tAugustusTrain\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={gid}.t1\n")
                out.write("###\n")

    write_gff3(out_train, train_set, 1)
    write_gff3(out_test, test_set, len(train_set) + 1)

    print(
        f"complete(>= {MIN_CDS}nt, repeat<{REPEAT_FRAC})={len(selected)} "
        f"after redundancy={len(kept)} train={len(train_set)} test={len(test_set)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
