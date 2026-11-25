#!/usr/bin/env bash
# PURPOSE: Build ONT BAM against review4 genome
# INPUT:   cho_v3_review4.fa + ${ONT_UL_1} ${ONT_UL_2}
# OUTPUT:  ont_to_review4.bam
# TOOL:    minimap2 map-ont + samtools
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -e
conda activate tgsgapcloser-env

REVIEW4="${WORK_DIR}/yahs/output_v3/juicer_post_review4/cho_v3_review4.fa"
ONT_1="${ONT_UL_1}"
ONT_2="${ONT_UL_2}"
OUTDIR="${WORK_DIR}/gapfill_review4"
BAM="${OUTDIR}/ont_to_review4.bam"
THREADS=16

mkdir -p "${OUTDIR}"

echo "===== $(date) Build ONT → review4 BAM for Flye gapfilling ====="
echo "Ref:    ${REVIEW4} ($(du -h ${REVIEW4} | cut -f1))"
echo "ONT:    ${ONT_1} ($(du -h ${ONT_1} | cut -f1))"
echo "        ${ONT_2} ($(du -h ${ONT_2} | cut -f1))"
echo "BAM:    ${BAM}"
echo "Threads: ${THREADS}"

echo ""
echo "[$(date)] Step 1: minimap2 + samtools sort..."
minimap2 -ax map-ont -t ${THREADS} "${REVIEW4}" "${ONT_1}" "${ONT_2}" \
    | samtools sort -@ 8 -o "${BAM}" -

echo "[$(date)] Step 2: samtools index..."
samtools index -@ 8 "${BAM}"

echo ""
echo "===== $(date) Done ====="
ls -lh "${BAM}" "${BAM}.bai" "${BAM}.csi" 2>/dev/null
