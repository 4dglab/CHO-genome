#!/usr/bin/env bash
# PURPOSE: Generate Juicebox .hic interaction file
# INPUT:   cho_yahs_v3_scaffolds_final.agp + hic_aligned.bam
# OUTPUT:  cho_yahs_v3.JBAT.hic
# TOOL:    juicer / juicer_tools 1.22.01
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -e

conda activate yahs-env

OUTDIR="${WORK_DIR}/yahs/output_v3"
BIN="${OUTDIR}/cho_yahs_v3.bin"
AGP="${OUTDIR}/cho_yahs_v3_scaffolds_final.agp"
FAI="${WORK_DIR}/purge_dups_v3/purged_clean.fa.fai"
JUICER_TOOLS="${JUICER_TOOLS}"
MAPQ=10
PREFIX="${OUTDIR}/cho_yahs_v3.JBAT"

cd "${OUTDIR}"

LOG="${OUTDIR}/jbhat_hic.log"
exec > >(tee -a "${LOG}") 2>&1

echo "===== $(date) Generate .JBAT.hic via SOP ====="
echo "BIN:    ${BIN} ($(du -h ${BIN} | cut -f1))"
echo "AGP:    ${AGP}"
echo "FAI:    ${FAI}"
echo "MAPQ:   ${MAPQ}"
echo "Prefix: ${PREFIX}"

# Step 1: juicer pre (input .bin, NOT bam!)
echo ""
echo "===== $(date) Step 1: juicer pre (from .bin) ====="

juicer pre \
    -a \
    -q "${MAPQ}" \
    -o "${PREFIX}" \
    "${BIN}" \
    "${AGP}" \
    "${FAI}" \
    > "${PREFIX}.log" \
    2>&1

echo "Step 1 done"
ls -lh "${PREFIX}.txt" "${PREFIX}.assembly" "${PREFIX}.assembly.agp" "${PREFIX}.liftover.agp" 2>/dev/null

# Step 2: Extract JBAT_SIZE from log
echo ""
echo "===== $(date) Step 2: Extract JBAT_SIZE ====="

JBAT_SIZE=$(awk '/^PRE_C_SIZE:[[:space:]]+assembly[[:space:]]+/ {print $3}' \
    "${PREFIX}.log" | tail -n 1)

if [ -z "${JBAT_SIZE}" ]; then
    echo "ERROR: Could not extract JBAT_SIZE from ${PREFIX}.log"
    echo "Check the log for PRE_C_SIZE line."
    grep "PRE_C_SIZE" "${PREFIX}.log" || true
    exit 1
fi

echo "JBAT_SIZE: assembly ${JBAT_SIZE}"

# Step 3: juicer_tools → .hic
echo ""
echo "===== $(date) Step 3: juicer_tools → .hic ====="

java -Xmx500g \
    -XX:ParallelGCThreads=4 \
    -jar "${JUICER_TOOLS}" \
    pre \
    --threads 8 \
    "${PREFIX}.txt" \
    "${PREFIX}.hic" \
    <(printf "assembly %s\n" "${JBAT_SIZE}") \
    > "${PREFIX}.hic.log" \
    2>&1

echo "Output: ${PREFIX}.hic ($(du -h ${PREFIX}.hic | cut -f1))"

echo ""
echo "===== $(date) Done ====="
ls -lh "${PREFIX}.txt" "${PREFIX}.assembly" "${PREFIX}.liftover.agp" "${PREFIX}.hic" 2>/dev/null
