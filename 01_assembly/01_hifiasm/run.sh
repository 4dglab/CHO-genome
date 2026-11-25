#!/usr/bin/env bash
# PURPOSE: hifiasm de novo assembly of HiFi reads (primary contigs)
# INPUT:   ${WORK_DIR}/hifi_fastq/CH6F_G1-{1,2}_hifi.fq.gz
# OUTPUT:  CH6F_G1_hifiasm.bp.p_ctg.fasta
# TOOL:    hifiasm v0.25.0
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# hifiasm HiFi genome assembly for Chinese hamster
# HiFi data: CH6F_G1-1_hifi.fq.gz, CH6F_G1-2_hifi.fq.gz

set -euo pipefail

# ---------- config ----------
THREADS=20
HIFI_DIR="${WORK_DIR}/hifi_fastq"
OUT_DIR="${WORK_DIR}/hifiasm_assembly"
OUT_PREFIX="${OUT_DIR}/CH6F_G1_hifiasm"

# ---------- env ----------
conda activate hifiasm-env

echo "============================================"
date
echo "Running hifiasm assembly..."
echo "Output prefix: ${OUT_PREFIX}"
echo "Threads: ${THREADS}"
echo "HiFi reads:"
echo "  ${HIFI_DIR}/CH6F_G1-1_hifi.fq.gz"
echo "  ${HIFI_DIR}/CH6F_G1-2_hifi.fq.gz"
echo "============================================"

hifiasm \
    -o "${OUT_PREFIX}" \
    -t "${THREADS}" \
    "${HIFI_DIR}/CH6F_G1-1_hifi.fq.gz" \
    "${HIFI_DIR}/CH6F_G1-2_hifi.fq.gz"

echo ""
echo "============================================"
date
echo "hifiasm assembly finished."
echo "Output files: ${OUT_DIR}"
echo "============================================"
