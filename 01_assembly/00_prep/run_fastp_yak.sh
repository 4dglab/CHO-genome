#!/usr/bin/env bash
# PURPOSE: Pre-process NGS and Hi-C reads (fastp QC/trimming) and build yak k-mer databases
# INPUT:   ${NGS_R1} ${NGS_R2} ${HIC_LIB1_R1} ${HIC_LIB1_R2} ${HIC_LIB2_R1} ${HIC_LIB2_R2}
# OUTPUT:  ${WORK_DIR}/nextpolish2/ngs_clean/R{1,2}.clean.fq.gz, ${WORK_DIR}/nextpolish2/kmer/k21.yak k31.yak,
#          ${WORK_DIR}/hi-c_data/clean/CH6F_G3-1-{1,2}_{R1,R2}.fq.gz
# TOOL:    fastp / yak
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/config.sh"

NGSCLEAN="${WORK_DIR}/nextpolish2/ngs_clean"
KMER="${WORK_DIR}/nextpolish2/kmer"
HICCLEAN="${WORK_DIR}/hi-c_data/clean"
mkdir -p "${NGSCLEAN}" "${KMER}" "${HICCLEAN}"

# ---- NGS: fastp QC/trimming ----
conda activate nextpolish2-env
fastp \
    -5 -3 -n 0 -f 5 -F 5 -t 5 -T 5 -q 20 \
    -i "${NGS_R1}" -I "${NGS_R2}" \
    -o "${NGSCLEAN}/R1.clean.fq.gz" -O "${NGSCLEAN}/R2.clean.fq.gz" \
    --thread "${THREADS}"

# ---- NGS: yak k-mer databases (reused by NextPolish2) ----
yak count -o "${KMER}/k21.yak" -k 21 -b 37 \
    <(zcat "${NGSCLEAN}/R1.clean.fq.gz") <(zcat "${NGSCLEAN}/R2.clean.fq.gz")
yak count -o "${KMER}/k31.yak" -k 31 -b 37 \
    <(zcat "${NGSCLEAN}/R1.clean.fq.gz") <(zcat "${NGSCLEAN}/R2.clean.fq.gz")

# ---- Hi-C: fastp QC/trimming (two libraries) ----
fastp -5 -3 -q 20 -i "${HIC_LIB1_R1}" -I "${HIC_LIB1_R2}" \
    -o "${HICCLEAN}/CH6F_G3-1-1_R1.fq.gz" -O "${HICCLEAN}/CH6F_G3-1-1_R2.fq.gz" \
    --thread "${THREADS}"
fastp -5 -3 -q 20 -i "${HIC_LIB2_R1}" -I "${HIC_LIB2_R2}" \
    -o "${HICCLEAN}/CH6F_G3-1-2_R1.fq.gz" -O "${HICCLEAN}/CH6F_G3-1-2_R2.fq.gz" \
    --thread "${THREADS}"

echo "Done: NGS clean + yak k-mer DBs + Hi-C clean reads"
