#!/usr/bin/env bash
# PURPOSE: Compute A/B compartments at 1 mb and 100 kb from the KR-balanced matrices.
# INPUT:   ${HI_C_DIR}/AB_compartments/K1_{1mb,100kb}_balanced.hic, ${GENOME_FA}
# OUTPUT:  ${HI_C_DIR}/AB_compartments/K1.{1mb,100kb}_balanced.{ab,bed,ev.txt}
# TOOL:    FAN-C v0.9.28 (fanc compartments; GC-corrected sign orientation)
# env: fanc-env
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   GENOME_FA          reference genome fasta with .fai index
#   FANC_BIN           fanc executable (fanc-env)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

ABDIR="${HI_C_DIR}/AB_compartments"
cd "${ABDIR}"

for RES in 1mb 100kb; do
    HIC="${ABDIR}/K1_${RES}_balanced.hic"
    echo "===== $(date) compartments ${RES} ====="
    "${FANC_BIN}" compartments \
        -d "K1.${RES}_balanced.bed" \
        -v "K1.${RES}_balanced.ev.txt" \
        -g "${GENOME_FA}" \
        "${HIC}" \
        "K1.${RES}_balanced.ab" \
        > "run_compartments_${RES}_balanced.log" 2>&1
    echo "DONE ${RES}"
done
echo "ALL_COMPARTMENTS_DONE"
ls -lh K1.*_balanced.bed
