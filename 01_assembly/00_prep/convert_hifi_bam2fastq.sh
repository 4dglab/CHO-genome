#!/usr/bin/env bash
# PURPOSE: Convert HiFi BAM to FASTQ.gz
# INPUT:   ${HIFI_BAM_1} ${HIFI_BAM_2}
# OUTPUT:  ${WORK_DIR}/hifi_fastq/CH6F_G1-{1,2}_hifi.fq.gz
# TOOL:    samtools fastq + pigz
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -e

conda activate verkko-env

THREADS=96
OUTDIR="${WORK_DIR}/hifi_fastq"

BAM_1="${HIFI_BAM_1}"
BAM_2="${HIFI_BAM_2}"

mkdir -p "${OUTDIR}"

# ====== BAM → FASTQ.gz ======
echo "===== $(date) Converting CH6F_G1-1 ====="
samtools fastq -T '*' "${BAM_1}" \
    | pigz -p ${THREADS} > "${OUTDIR}/CH6F_G1-1_hifi.fq.gz"
echo "===== $(date) CH6F_G1-1 done ====="

echo "===== $(date) Converting CH6F_G1-2 ====="
samtools fastq -T '*' "${BAM_2}" \
    | pigz -p ${THREADS} > "${OUTDIR}/CH6F_G1-2_hifi.fq.gz"
echo "===== $(date) CH6F_G1-2 done ====="

echo "===== $(date) All conversions finished ====="
