#!/usr/bin/env bash
# PURPOSE: Flye gap-filling (calls fill_gaps_flye.py)
# INPUT:   cho_v3_review4.fa + ont_to_review4.bam
# OUTPUT:  cho_final_flye_gapfilled.fa
# TOOL:    Flye 2.9.6
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -e
conda activate gapfill-env

REF="${WORK_DIR}/yahs/output_v3/juicer_post_review4/cho_v3_review4.fa"
BAM="${WORK_DIR}/gapfill_review4/ont_to_review4.bam"
OUTDIR="${WORK_DIR}/gapfill_review4/flye_output"
THREADS=4

echo "===== $(date) Flye gapfilling RETRY (TIMEOUT=60min, MIN_READS=2) ====="
echo "Ref:     ${REF}  ($(du -h "${REF}" | cut -f1))"
echo "BAM:     ${BAM}  ($(du -h "${BAM}" | cut -f1))"
echo "Outdir:  ${OUTDIR}"
echo "Workers: ${THREADS}"
echo ""

python3 ${WORK_DIR}/gapfill_review4/fill_gaps_flye_retry.py \
    -r "${REF}" \
    -b "${BAM}" \
    -o "${OUTDIR}" \
    -p ${THREADS} \
    --resume

echo ""
echo "===== $(date) Done ====="
ls -lh "${OUTDIR}/cho_final_flye_gapfilled.fa" 2>/dev/null
