#!/usr/bin/env python3
# PURPOSE: Flye per-gap filling core logic
# INPUT:   reference genome + ONT BAM
# OUTPUT:  cho_final_flye_gapfilled.fa
# TOOL:    Flye / minimap2 / samtools


"""
Batch gap-filling using Flye de novo assembly of ONT reads.
Method:
  1. Extract flank reference (with N-gaps) around each gap.
  2. Pull ONT reads spanning the gap from pre-aligned BAM.
  3. Run Flye (--min-overlap 1000 --read-error 0.1) on those reads.
  4. Align Flye contigs to the flank reference (with N's) via minimap2.
  5. Parse CIGAR to locate query bases mapping to the reference N-region.
  6. Extract fill sequence and replace N's in the genome.

References:
  - Flye: Kolmogorov et al., 2019, Nature Biotechnology
  - minimap2: Li, 2018, Bioinformatics
"""

import subprocess, os, sys, time, argparse, re
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

# ── absolute binary paths (no conda env required in subprocess) ──────────
CONDA_PREFIX = os.environ.get("CONDA_PREFIX_GAPFILL", "/path/to/anaconda3/envs/gapfill-env")
SAMTOOLS   = f"{CONDA_PREFIX}/bin/samtools"
MINIMAP2   = f"{CONDA_PREFIX}/bin/minimap2"
FLYE       = f"{CONDA_PREFIX}/bin/flye"

FLYE_THREADS = 4
FLYE_TIMEOUT = 3600      # 60 min per gap (retry)
MIN_READS    = 2          # lowered for retry
SMALL_GAP    = 10        # skip tiny N runs (< 10 bp)

# ── helpers ────────────────────────────────────────────────────────────────

def load_fasta(path):
    seqs, cur = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, name, seq, width=60):
    """Write a single-entry FASTA file."""
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i : i + width] + "\n")


def run_cmd(cmd, timeout=120, check=False):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=timeout,
        env={**os.environ, "PATH": f"{CONDA_PREFIX}/bin:{os.environ.get('PATH', '')}"}
    )
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {cmd}\n{result.stderr[:500]}")
    return result.stdout.strip(), result.stderr.strip(), result.returncode


def count_fastq(path):
    try:
        with open(path) as f:
            return sum(1 for _ in f) // 4
    except FileNotFoundError:
        return 0


# ── per-gap logic ──────────────────────────────────────────────────────────

def fill_one_gap(args):
    """
    Process one gap in a subprocess worker.
    args = (scaffold, gs, gap_len, ref_path, bam_path, flank, work_dir, scaffold_len)

    Returns (gap_id, success: bool, message: str, fill_seq: str or None).
    """
    scaffold, gs, gap_len, ref_path, bam_path, flank, work_dir, scaffold_len = args
    gap_id = f"{scaffold}_{gs}"

    wdir = Path(work_dir) / gap_id
    wdir.mkdir(parents=True, exist_ok=True)

    ge_incl  = gs + gap_len - 1     # 1-based inclusive end
    ge_excl  = gs + gap_len          # 1-based exclusive end (= first bp after gap)
    left     = max(gs - flank, 1)
    right    = min(ge_incl + flank, scaffold_len)
    region   = f"{scaffold}:{left}-{right}"

    try:
        # ── 1. extract flank reference (keep N's so we can locate the gap) ──
        _, _, rc = run_cmd(
            f"{SAMTOOLS} faidx '{ref_path}' '{region}' 2>/dev/null | "
            f"grep -v '^>' | tr -d '\\n' > {wdir}/flank_ref.txt",
            timeout=60,
        )
        with open(wdir / "flank_ref.txt") as fh:
            flank_seq = fh.read().strip()
        if not flank_seq or rc != 0:
            return (gap_id, False, "empty flank reference", None)

        # verify gap position inside the extracted flank
        # flank_seq extracted from genome position 'left' (1-based)
        # gap starts at genome position 'gs' (1-based)
        # gap offset in flank_seq = gs - left
        ref_gap_start = gs - left          # 0-based offset inside flank_seq
        ref_gap_end   = ref_gap_start + gap_len
        if flank_seq[ref_gap_start:ref_gap_end] != "N" * gap_len:
            return (gap_id, False, "gap N positions unexpected", None)

        # write flank reference as FASTA (with N's)
        write_fasta(wdir / "flank_ref.fa", region, flank_seq)

        # ── 2. extract cross-gap ONT reads ──
        # If gap_reads.fastq already exists (extracted by awk from minimap2 SAM),
        # skip extraction. Otherwise, fall back to BAM-based extraction with
        # CIGAR-corrected reference span computation.
        gap_reads = wdir / "gap_reads.fastq"
        if gap_reads.exists() and gap_reads.stat().st_size > 0:
            pass  # reads already extracted by upstream awk script
        else:
            awk_file = wdir / "filter.awk"
            awk_content = (
                "/^@/ { print; next }\n"
                "{\n"
                "    s = 0\n"
                "    c = $6\n"
                "    while (match(c, /[0-9]+[MIDNSHP=X]/)) {\n"
                "        m = substr(c, RSTART, RLENGTH)\n"
                "        l = m + 0\n"
                "        op = substr(m, RLENGTH, 1)\n"
                "        if (op == \"M\" || op == \"D\" || op == \"N\" || op == \"=\" || op == \"X\") {\n"
                "            s += l\n"
                "        }\n"
                "        c = substr(c, RSTART + RLENGTH)\n"
                "    }\n"
                f"    re = $4 + s - 1\n"
                f"    if (and($2, 4) == 0 && $4 <= {gs} && re >= {ge_incl}) {{\n"
                "        print\n"
                "    }\n"
                "}\n"
            )
            with open(awk_file, "w") as af:
                af.write(awk_content)
            _, _, _ = run_cmd(
                f"{SAMTOOLS} view -h '{bam_path}' '{region}' 2>/dev/null | "
                f"awk -f {awk_file} | "
                f"{SAMTOOLS} fastq - 2>/dev/null > {gap_reads}",
                timeout=180,
            )

        n_reads = count_fastq(gap_reads)
        if n_reads < MIN_READS:
            return (gap_id, False, f"insufficient reads ({n_reads})", None)

# Deduplicate and subsample to MAX_READS
        MAX_READS = 500
        flye_reads = gap_reads
        if n_reads > MAX_READS:
            flye_reads = wdir / "gap_reads_subsampled.fastq"
            if not flye_reads.exists() or count_fastq(flye_reads) < MIN_READS:
                # Deduplicate by read name, then take first MAX_READS
                seen = set()
                count = 0
                with open(gap_reads) as fin, open(flye_reads, "w") as fout:
                    while count < MAX_READS:
                        h = fin.readline()
                        if not h:
                            break
                        s = fin.readline()
                        p = fin.readline()
                        q = fin.readline()
                        if not q:
                            break
                        rid = h.strip().split()[0]  # @readname
                        if rid not in seen:
                            seen.add(rid)
                            fout.write(h + s + p + q)
                            count += 1
            sub_n = count_fastq(flye_reads)
            if sub_n < MIN_READS:
                return (gap_id, False, f"subsampled to {sub_n} reads", None)
        else:
            # Deduplicate even for small read sets (Flye requires unique IDs)
            dedup = wdir / "gap_reads_dedup.fastq"
            if not dedup.exists() or count_fastq(dedup) != n_reads:
                seen = set()
                count = 0
                with open(gap_reads) as fin, open(dedup, "w") as fout:
                    while True:
                        h = fin.readline()
                        if not h:
                            break
                        s = fin.readline()
                        p = fin.readline()
                        q = fin.readline()
                        if not q:
                            break
                        rid = h.strip().split()[0]
                        if rid not in seen:
                            seen.add(rid)
                            fout.write(h + s + p + q)
                            count += 1
                dedup_n = count_fastq(dedup)
                if dedup_n < MIN_READS:
                    dedup.unlink(missing_ok=True)
                    return (gap_id, False, f"after dedup only {dedup_n} reads", None)
                flye_reads = dedup

        # ── 3. Flye de-novo assembly ──
        flye_dir = wdir / "flye_out"
        flye_dir.mkdir(exist_ok=True)
        flye_cmd = (
            f"{FLYE} --nano-hq {flye_reads} "
            f"-o {flye_dir} -g 100k --min-overlap 1000 "
            f"--read-error 0.1 -t {FLYE_THREADS} --asm-coverage 40"
        )
        stdout, stderr, rc = run_cmd(flye_cmd, timeout=FLYE_TIMEOUT)
        if rc != 0:
            return (gap_id, False, f"Flye failed (rc={rc})", None)

        flye_asm = flye_dir / "assembly.fasta"
        if not flye_asm.exists():
            return (gap_id, False, "no Flye assembly.fasta", None)

        # ── 4. align Flye contigs to flank reference (with N's) ──
        paf_out, _, _ = run_cmd(
            f"{MINIMAP2} -t 4 -c {wdir}/flank_ref.fa {flye_asm} 2>/dev/null",
            timeout=120,
        )
        if not paf_out:
            return (gap_id, False, "no Flye contig alignment to flank", None)

        # ── 5. find best contig that spans the N-gap region ──
        best_contig  = None
        best_score   = 0.0
        best_paf     = None
        flank_len    = len(flank_seq)
        ref_target   = ref_gap_start           # gap starts here in flank ref (0-based)
        ref_gap_span = ref_gap_end + 1         # one past gap end

        for line in paf_out.split("\n"):
            line = line.strip()
            if not line:
                continue
            flds = line.split("\t")
            if len(flds) < 12:
                continue
            qname  = flds[0]
            qlen   = int(flds[1])
            qstart = int(flds[2])
            qend   = int(flds[3])
            tname  = flds[5]
            tlen   = int(flds[6])
            tstart = int(flds[7])
            tend   = int(flds[8])
            matches    = int(flds[9])
            aln_block  = int(flds[10])

            # check that alignment spans the N region
            if not (tstart < ref_target and tend > ref_gap_end):
                continue

            identity = matches / max(aln_block, 1)
            span     = tend - tstart
            score    = span * identity
            if score > best_score:
                best_score   = score
                best_contig  = qname
                best_paf     = flds

        if best_contig is None:
            return (gap_id, False, "no contig spans the gap", None)

        # ── 6. extract gap-fill sequence via CIGAR ──
        flye_contigs = load_fasta(flye_asm)
        contig_seq   = flye_contigs.get(best_contig, "")
        if not contig_seq:
            return (gap_id, False, f"contig {best_contig} not found in FASTA", None)

        qstart = int(best_paf[2])
        tstart = int(best_paf[7])

        # try cg:Z (CIGAR) tag
        cg_tag = ""
        for tag in best_paf[11:]:
            if tag.startswith("cg:Z:"):
                cg_tag = tag[5:]
                break

        fill_seq = ""
        if cg_tag:
            # parse CIGAR to map ref N-positions → query positions
            ops = re.findall(r"(\d+)([MID])", cg_tag)
            ref_pos  = tstart
            query_pos = qstart
            gap_query_start = None
            gap_query_end   = None

            for length_str, op in ops:
                length = int(length_str)
                if op == "M":
                    if ref_pos <= ref_target < ref_pos + length:
                        gap_query_start = query_pos + (ref_target - ref_pos)
                    if ref_pos <= ref_gap_end <= ref_pos + length:
                        gap_query_end = query_pos + (ref_gap_end - ref_pos)
                        break
                    ref_pos   += length
                    query_pos += length
                elif op == "I":
                    query_pos += length
                elif op == "D":
                    ref_pos   += length
        else:
            # no CIGAR – approximate: gap N region is roughly in the middle of the alignment
            aligned_region = contig_seq[qstart:int(best_paf[3])]
            ratio_start = (ref_target - tstart) / max(int(best_paf[8]) - tstart, 1)
            ratio_end   = (ref_gap_end  - tstart) / max(int(best_paf[8]) - tstart, 1)
            gap_query_start = int(len(aligned_region) * ratio_start)
            gap_query_end   = int(len(aligned_region) * ratio_end)

        if gap_query_start is None or gap_query_end is None:
            return (gap_id, False, "could not map gap in CIGAR", None)

        raw_fill = contig_seq[gap_query_start:gap_query_end]
        fill_n   = raw_fill.count("N")

        # quality check: at least half must be non-N
        if not raw_fill or fill_n >= len(raw_fill) * 0.5:
            return (gap_id, False, f"bad fill ({len(raw_fill)} bp, {100*fill_n/max(len(raw_fill),1):.0f}% N)", None)

        # fit to expected gap length
        if len(raw_fill) < gap_len:
            fill_seq = raw_fill + "N" * (gap_len - len(raw_fill))
        else:
            # take central gap_len bases
            excess = len(raw_fill) - gap_len
            start  = excess // 2
            fill_seq = raw_fill[start : start + gap_len]

        # final sanity
        if fill_seq.count("N") > gap_len * 0.5:
            return (gap_id, False, f"still too many Ns after fitting", None)

        with open(wdir / "fill_seq.txt", "w") as fh:
            fh.write(fill_seq)

        return (gap_id, True,
                f"OK reads={n_reads} contig={best_contig} fill={len(raw_fill)}→{len(fill_seq)}bp N={fill_n}",
                fill_seq)

    except subprocess.TimeoutExpired:
        return (gap_id, False, "timeout", None)
    except Exception:
        return (gap_id, False, traceback.format_exc()[-200:], None)


# ── main ───────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Batch gap-filling with Flye + ONT reads")
    ap.add_argument("-r", "--ref",    required=True, help="Reference FASTA with N-gaps")
    ap.add_argument("-b", "--bam",    required=True, help="ONT reads aligned to reference (BAM)")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument("-f", "--flank",  type=int, default=2000, help="Flanking bp each side [2000]")
    ap.add_argument("-p", "--parallel", type=int, default=4, help="Parallel workers [4]")
    ap.add_argument("--resume",       action="store_true", help="Skip gaps with fill_seq.txt already present")
    ap.add_argument("--dry-run",      action="store_true", help="Only report gaps, exit")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] Loading genome from {args.ref}")
    genome = load_fasta(args.ref)
    scaffold_lengths = {s: len(seq) for s, seq in genome.items()}
    total_N_before = sum(seq.count("N") for seq in genome.values())
    print(f"  scaffolds={len(genome)}  total_N={total_N_before}")

    # discover gaps
    gaps = []
    for scaffold, seq in sorted(genome.items(), key=lambda x: len(x[1]), reverse=True):
        for m in re.finditer(r"N+", seq):
            glen = m.end() - m.start()
            if glen < SMALL_GAP:
                continue
            gs = m.start() + 1            # 1-based
            gaps.append((scaffold, gs, glen))
    print(f"  gaps to fill = {len(gaps)}")

    if args.dry_run:
        return

    # resume logic
    if args.resume:
        todo = []
        skipped = 0
        for g in gaps:
            d = outdir / f"{g[0]}_{g[1]}"
            fpath = d / "fill_seq.txt"
            if fpath.exists():
                seq = fpath.read_text().strip()
                if seq and seq.count("N") < len(seq) * 0.5:
                    s, p, l = g
                    genome[s] = genome[s][:p-1] + seq + genome[s][p-1+l:]
                    skipped += 1
                    continue
            todo.append(g)
        print(f"  resumed: {skipped} already done, {len(todo)} remaining")
        gaps = todo

    if not gaps:
        print("Nothing to do.")
        return

    # prepare worker arguments
    worker_args = [
        (scaffold, gs, glen, args.ref, args.bam, args.flank, str(outdir), scaffold_lengths[scaffold])
        for scaffold, gs, glen in gaps
    ]

    t0 = time.time()
    filled, failed = 0, 0

    print(f"\n[{time.strftime('%H:%M:%S')}] Processing {len(gaps)} gaps with {args.parallel} workers...\n")
    sys.stdout.flush()

    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(fill_one_gap, a): a for a in worker_args}
        for fut in as_completed(futures):
            gap_id, ok, msg, fill_seq = fut.result()

            # map back to gap info
            scaffold = gap_id.rsplit("_", 1)[0]
            gs_str   = gap_id.rsplit("_", 1)[1]
            gs       = int(gs_str)
            gap_len  = None
            for s, p, l in gaps:
                if s == scaffold and p == gs:
                    gap_len = l
                    break

            if ok:
                genome[scaffold] = (
                    genome[scaffold][:gs-1] + fill_seq + genome[scaffold][gs-1+gap_len:]
                )
                filled += 1
            else:
                failed += 1

            elapsed = time.time() - t0
            eta = elapsed / max(filled + failed, 1) * (len(gaps) - filled - failed)
            print(f"[{time.strftime('%H:%M:%S')}] {gap_id:>25s}  {msg:70s}  "
                  f"f={filled} fail={failed}  {elapsed:.0f}s  eta={eta/60:.0f}min")
            sys.stdout.flush()

    # write final genome
    output_fa = outdir / "cho_final_flye_gapfilled.fa"
    print(f"\n[{time.strftime('%H:%M:%S')}] Writing {output_fa} ...")
    with open(output_fa, "w") as fh:
        for scaffold in sorted(genome.keys(), key=lambda s: scaffold_lengths[s], reverse=True):
            seq = genome[scaffold]
            fh.write(f">{scaffold}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    total_N_after = sum(seq.count("N") for seq in genome.values())
    elapsed = time.time() - t0

    print(f"""
{'='*60}
FLYE GAP-FILLING REPORT
{'='*60}
Gaps:       {len(gaps)}
Filled:     {filled}  ({100*filled/len(gaps):.1f}%)
Failed:     {failed}
N before:   {total_N_before} bp
N after:    {total_N_after} bp
N removed:  {total_N_before - total_N_after} bp
Time:       {elapsed:.0f}s ({elapsed/3600:.1f}h)
Output:     {output_fa}
{'='*60}
""")


if __name__ == "__main__":
    main()
