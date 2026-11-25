#!/usr/bin/env python3
# PURPOSE: Merge HiTE+Dfam repeat GFF3
# INPUT:   HiTE/Dfam GFF3
# OUTPUT:  merge_all.repeat.gff3
# TOOL:    custom


"""Merge HiTE de novo + Dfam RepeatMasker repeat annotations into one RepeatMasker-style GFF3.

Reads the native RepeatMasker .out tables (HiTE.out + genome.fa.out) and converts
each entry into a GFF line in RepeatMasker .out.gff format:

    <seqid>  RepeatMasker  <root_class>  <qbegin>  <qend>  <SW>  <strand>  .  ID=TE%07d;Target=<family> <rbegin> <rend>;Class=<root_class>;PercDiv=<div>;PercDel=<del>;PercIns=<ins>;

Rules:
- feature type / Class = root of the repeat class/family (text after last '/'; e.g. LINE/L1 -> L1).
- ALL repeats are kept (Simple_repeat / Low_complexity / satellite included).
- Dedup: a Dfam feature is dropped if >= 90% of its length is covered by a HiTE feature
  (de novo full-length TE spans take priority).
- IDs are re-numbered sequentially TE0000001..N across the merged, position-sorted output.
- A ##sequence-region header line is emitted for EVERY contig in the genome (multi-contig aware).
"""
import argparse
import re
import sys


def parse_out(fpath):
    feats = []
    with open(fpath) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip() or line.startswith(' '):
                if line.lstrip().startswith('SW') or line.lstrip().startswith('score'):
                    continue
            tokens = [t for t in line.split() if t]
            if len(tokens) < 15:
                continue
            if tokens[-1] == '*':
                tokens = tokens[:-1]
            if len(tokens) != 15:
                sys.stderr.write('skip malformed line: {}\n'.format(line))
                continue

            def coord(x):
                if x.startswith('(') and x.endswith(')'):
                    return -int(x[1:-1])
                return int(x)

            feats.append({
                'sw': tokens[0], 'div': tokens[1], 'del': tokens[2], 'ins': tokens[3],
                'seqid': tokens[4],
                'qbegin': int(tokens[5]), 'qend': int(tokens[6]),
                'strand': tokens[8],
                'family': tokens[9], 'class': tokens[10],
                'rbegin': coord(tokens[11]), 'rend': coord(tokens[12]),
            })
    return feats


def genome_lengths(fa_path):
    lengths = {}
    seqid = None
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith('>'):
                seqid = line.split()[0][1:]
                lengths[seqid] = 0
            elif seqid is not None:
                lengths[seqid] += len(line.strip())
    return lengths


def root_class(cls):
    return cls.split('/')[-1] if '/' in cls else cls


def format_line(f, idx):
    root = root_class(f['class'])
    strand = '-' if f['strand'] == 'C' else f['strand']
    attr = 'ID=TE{:07d};Target={} {} {};Class={};PercDiv={};PercDel={};PercIns={};'.format(
        idx, f['family'], f['rbegin'], f['rend'], f['class'], f['div'], f['del'], f['ins'])
    return '\t'.join([
        f['seqid'], 'RepeatMasker', root,
        str(f['qbegin']), str(f['qend']), f['sw'], strand, '.', attr,
    ])


def cov_frac(f, h):
    s = max(f['qbegin'], h['qbegin'])
    e = min(f['qend'], h['qend'])
    if e < s:
        return 0.0
    overlap = e - s + 1
    return overlap / (f['qend'] - f['qbegin'] + 1)


def dedup_dfam_vs_hite(dfam, hite):
    hite = sorted(hite, key=lambda x: (x['seqid'], x['qbegin'], x['qend']))
    dfam_sorted = sorted(dfam, key=lambda x: (x['seqid'], x['qbegin'], x['qend']))

    kept = []
    dropped = 0
    pool = []
    hidx = 0
    for f in dfam_sorted:
        while hidx < len(hite) and hite[hidx]['seqid'] == f['seqid'] and hite[hidx]['qbegin'] <= f['qend']:
            pool.append(hite[hidx])
            hidx += 1
        pool = [h for h in pool if h['qend'] >= f['qbegin']]
        if any(cov_frac(f, h) >= 0.9 for h in pool):
            dropped += 1
            continue
        kept.append(f)
    return kept, dropped


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hite', required=True, help='HiTE .out table')
    ap.add_argument('--dfam', required=True, help='Dfam RepeatMasker .out table')
    ap.add_argument('--genome', required=True, help='genome FASTA used to emit per-contig ##sequence-region headers')
    ap.add_argument('--out', required=True, help='output merge_all.repeat.gff3')
    args = ap.parse_args()

    hite = parse_out(args.hite)
    dfam = parse_out(args.dfam)

    dfam_kept, dropped = dedup_dfam_vs_hite(dfam, hite)

    lines = [format_line(f, i + 1) for i, f in enumerate(hite)]
    base = len(hite)
    lines += [format_line(f, base + i + 1) for i, f in enumerate(dfam_kept)]

    def sort_key(line):
        cols = line.split('\t')
        return (cols[0], int(cols[3]), int(cols[4]))
    lines.sort(key=sort_key)

    lengths = genome_lengths(args.genome)
    with open(args.out, 'w') as fh:
        fh.write('##gff-version 3\n')
        for seqid in lengths:
            fh.write('##sequence-region {} 1 {}\n'.format(seqid, lengths[seqid]))
        for line in lines:
            fh.write(line + '\n')

    print('=== merge_all.repeat.gff3 stats (RepeatMasker-style) ===')
    print('HiTE features      :', len(hite))
    print('Dfam features in   :', len(dfam))
    print('Dfam dropped (>=90% covered by HiTE):', dropped)
    print('Dfam kept          :', len(dfam_kept))
    print('total merged       :', len(lines))
    print('sequence-region    : {} contigs'.format(len(lengths)))
    print('written to         :', args.out)


if __name__ == '__main__':
    main()
