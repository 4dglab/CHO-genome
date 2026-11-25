#!/usr/bin/env python3
# PURPOSE: StringTie GTF -> GFF3 (multi-contig fix)
# INPUT:   transcripts_merged.gtf
# OUTPUT:  transcripts_merged.gff3
# TOOL:    custom


"""Convert StringTie GTF to GFF3 for Augustus training input.

Multi-contig aware (74-contig full genome).
"""
import sys
from collections import OrderedDict

def convert(gtf_path, gff3_path):
    transcripts = OrderedDict()   # transcript_id -> dict
    genes = OrderedDict()         # gene_id -> [transcript_ids]

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0]
            feat = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]

            gene_id = None
            transcript_id = None
            for a in attrs.rstrip(";").split("; "):
                if a.startswith("gene_id"):
                    gene_id = a.split('"')[1]
                elif a.startswith("transcript_id"):
                    transcript_id = a.split('"')[1]

            if feat == "transcript":
                transcripts[transcript_id] = {
                    "seqid": chrom, "start": start, "end": end, "strand": strand,
                    "gene_id": gene_id, "exons": []
                }
                if gene_id not in genes:
                    genes[gene_id] = []
                genes[gene_id].append(transcript_id)
            elif feat == "exon" and transcript_id and transcript_id in transcripts:
                transcripts[transcript_id]["exons"].append((start, end))

    with open(gff3_path, "w") as out:
        out.write("##gff-version 3\n")
        for gene_id, tids in genes.items():
            t_first = transcripts[tids[0]]
            g_seqid = t_first["seqid"]
            g_strand = t_first["strand"]
            g_starts = [transcripts[t]["start"] for t in tids]
            g_ends = [transcripts[t]["end"] for t in tids]
            out.write(f"{g_seqid}\tStringTie\tgene\t{min(g_starts)}\t{max(g_ends)}\t.\t{g_strand}\t.\tID={gene_id}\n")
            for tid in tids:
                t = transcripts[tid]
                out.write(f"{t['seqid']}\tStringTie\tmRNA\t{t['start']}\t{t['end']}\t.\t{t['strand']}\t.\tID={tid};Parent={gene_id}\n")
                for i, (es, ee) in enumerate(t["exons"], 1):
                    out.write(f"{t['seqid']}\tStringTie\texon\t{es}\t{ee}\t.\t{t['strand']}\t.\tID={tid}.exon{i};Parent={tid}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 gtf2gff3.py input.gtf output.gff3", file=sys.stderr)
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
    print(f"Converted {sys.argv[1]} -> {sys.argv[2]}")
