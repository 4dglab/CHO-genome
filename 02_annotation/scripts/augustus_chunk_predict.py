#!/usr/bin/env python3
# PURPOSE: Augustus chunked parallel ab initio prediction
# INPUT:   softmasked genome + training set
# OUTPUT:  augustus_abinitio.gff3
# TOOL:    Augustus 3.5.0


"""Run Augustus in parallel by splitting the genome into chunks.

Handles coordinate remapping and deduplication of genes that span chunk
boundaries (with overlap, each such gene is predicted in >=1 chunk; duplicate
copies with reciprocal CDS overlap >= 0.5 are merged keeping the longer gene).

Usage:
    augustus_chunk_predict.py genome.fa species out.gff3
        [--threads 16] [--chunk 5000000] [--overlap 200000]
        [--augustus /path/to/augustus]
"""

import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def split_fasta(path, chunk_size, overlap):
    """Yield (chunk_id, seqid, offset0based, start1based, end1based, seqstr)."""
    seqs = {}
    with open(path) as fh:
        cur = None
        buf = []
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line.split()[0][1:]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf)
    chunks = []
    for seqid, seq in seqs.items():
        L = len(seq)
        pos = 0
        idx = 0
        while pos < L:
            end = min(L, pos + chunk_size)
            chunk_start1 = pos + 1
            chunk_end1 = end
            chunks.append((f"{seqid}|chunk{idx}", seqid, pos, chunk_start1, chunk_end1, seq[pos:end]))
            pos += chunk_size - overlap
            idx += 1
    return chunks


def run_augustus(seqid, offset, seq, species, augustus, workdir):
    fa = os.path.join(workdir, f"{seqid.replace('|', '__')}.fa")
    with open(fa, "w") as fh:
        fh.write(f">{seqid}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i : i + 80] + "\n")
    cmd = [
        augustus, f"--species={species}", "--gff3=on", "--softmasking=1",
        "--UTR=off", "--print_utr=off", "--progress=false", fa,
    ]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as e:
        print(f"WARN: augustus failed for {seqid}: {e}", file=sys.stderr)
        print(e.stderr[-2000:], file=sys.stderr)
        return []
    genes = parse_augustus_gff(out)
    for g in genes:
        g["seqid"] = seqid
        g["offset"] = offset
    return genes


def parse_augustus_gff(text):
    genes = []
    cur_gene = None
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        ftype = cols[2]
        if ftype == "gene":
            cur_gene = {
                "seqid": cols[0],
                "type": "gene",
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "score": cols[5],
                "attrs": cols[8],
                "children": [],
            }
            genes.append(cur_gene)
        elif ftype in ("transcript", "CDS", "exon", "stop_codon", "start_codon", "intron", "UTR"):
            if cur_gene is not None:
                cur_gene["children"].append(
                    {
                        "type": ftype,
                        "start": int(cols[3]),
                        "end": int(cols[4]),
                        "strand": cols[6],
                        "phase": cols[7],
                        "attrs": cols[8],
                        "seqid": cols[0],
                    }
                )
    return genes


def remap(g, chunk_seqid):
    off = g["offset"]
    g["start"] += off
    g["end"] += off
    g["seqid"] = chunk_seqid
    for c in g["children"]:
        c["start"] += off
        c["end"] += off
        c["seqid"] = chunk_seqid
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_fa")
    ap.add_argument("species")
    ap.add_argument("out_gff")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--chunk", type=int, default=5000000)
    ap.add_argument("--overlap", type=int, default=200000)
    ap.add_argument("--augustus", default="augustus")
    args = ap.parse_args()

    chunks = split_fasta(args.genome_fa, args.chunk, args.overlap)

    all_genes = []
    with tempfile.TemporaryDirectory(prefix="augchunk_") as workdir:
        jobs = [
            (cid, seqid, off, start1, end1, seq)
            for (cid, seqid, off, start1, end1, seq) in chunks
        ]
        results = [None] * len(jobs)
        import concurrent.futures as cf

        def run(job):
            cid, seqid, off, start1, end1, seq = job
            return off, run_augustus(seqid, off, seq, args.species, args.augustus, workdir)

        with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(run, j) for j in jobs]
            for f in cf.as_completed(futs):
                off, gs = f.result()
                all_genes.extend(gs)

    for g in all_genes:
        remap(g, g["seqid"].split("|")[0])

    all_genes.sort(key=lambda g: (g["seqid"], g["start"]))

    # dedup overlapping genes (reciprocal >= 0.5), keep longer
    kept = []
    for g in all_genes:
        g_len = g["end"] - g["start"] + 1
        dup = False
        removals = []
        for i, h in enumerate(kept):
            if g["seqid"] != h["seqid"] or g["strand"] != h["strand"]:
                continue
            h_len = h["end"] - h["start"] + 1
            inter = max(0, min(g["end"], h["end"]) - max(g["start"], h["start"]) + 1)
            if inter / g_len >= 0.5 and inter / h_len >= 0.5:
                if g_len <= h_len:
                    dup = True
                    break
                removals.append(i)
        if dup:
            continue
        kept = [k for i, k in enumerate(kept) if i not in removals]
        kept.append(g)

    with open(args.out_gff, "w") as out:
        out.write("##gff-version 3\n")
        for gi, g in enumerate(kept, 1):
            gid = f"aug_g{gi}"
            out.write(
                f"{g['seqid']}\taugustus\tgene\t{g['start']}\t{g['end']}\t{g['score']}\t{g['strand']}\t.\tID={gid}\n"
            )
            transcripts = [c for c in g["children"] if c["type"] == "transcript"]
            cdss = [c for c in g["children"] if c["type"] == "CDS"]
            for ntr, c in enumerate(transcripts, 1):
                tid = f"{gid}.t{ntr}"
                out.write(
                    f"{c['seqid']}\taugustus\ttranscript\t{c['start']}\t{c['end']}\t{g['score']}\t{c['strand']}\t.\tID={tid};Parent={gid}\n"
                )
                c["tid"] = tid
            for c in sorted(cdss, key=lambda x: x["start"]):
                tid = c.get("tid") or (transcripts[0]["tid"] if transcripts else gid)
                out.write(
                    f"{c['seqid']}\taugustus\tCDS\t{c['start']}\t{c['end']}\t.\t{c['strand']}\t{c['phase']}\tParent={tid}\n"
                )
            out.write("###\n")

    n_genes = len(kept)
    print(f"predicted {n_genes} genes -> {args.out_gff}", file=sys.stderr)


if __name__ == "__main__":
    main()
