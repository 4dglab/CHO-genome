#!/usr/bin/env bash
# PURPOSE: Call chromatin loops with Juicer Tools HiCCUPS on the merged Hi-C matrix.
# INPUT:   ${HIC_FILE} (merged .hic)
# OUTPUT:  ${HI_C_DIR}/loop/K1_merged.hiccups/merged_loops.bedpe
# TOOL:    Juicer Tools 1.22.01 (hiccups, KR-normalized)
# env: hic-tools
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   HIC_FILE           merged .hic matrix (HiC-Pro valid pairs -> juicer pre)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

OUTDIR="${HI_C_DIR}/loop/K1_merged.hiccups"
mkdir -p "${OUTDIR}"
LOG="${HI_C_DIR}/loop/run_hiccups.log"
mkdir -p "$(dirname "${LOG}")"

exec > >(tee -a "${LOG}") 2>&1

echo "===== $(date) HiCCUPS loop calling ====="
echo "Input : ${HIC_FILE}"
echo "Output: ${OUTDIR}"
echo "Threads: ${THREADS}"

java -jar "${JUICER_TOOLS}" hiccups \
  --threads "${THREADS}" \
  -r 5000,10000 \
  "${HIC_FILE}" \
  "${OUTDIR}" \
  -k KR -m 500 -f 0.15,0.15 -p 6,3 -i 10,6 \
  -t 2,0.01,0.01,0.01 -d 10000,10000 --cpu

echo "===== $(date) done ====="
ls -lh "${OUTDIR}/merged_loops.bedpe"
