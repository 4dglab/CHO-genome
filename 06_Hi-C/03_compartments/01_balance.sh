#!/usr/bin/env bash
# PURPOSE: KR-balance the 1 mb and 100 kb Hi-C matrices and export balanced .hic files.
# INPUT:   ${MCCOOL} (@1mb and @100kb resolutions)
# OUTPUT:  ${HI_C_DIR}/AB_compartments/K1_{1mb,100kb}_balanced.hic
# TOOL:    FAN-C v0.9.28 (fanc hic -m KR)
# env: fanc-env
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   MCCOOL             merged .mcool matrix (hic2cool output)
#   FANC_BIN           fanc executable (fanc-env)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

ABDIR="${HI_C_DIR}/AB_compartments"
mkdir -p "${ABDIR}"
cd "${ABDIR}"

for RES in 1mb 100kb; do
    OUT="K1_${RES}_balanced.hic"
    if [ -s "${OUT}" ]; then
        echo "SKIP (exists): ${OUT}"
        continue
    fi
    echo "===== $(date) balance ${RES} ====="
    "${FANC_BIN}" hic \
        -b "${RES}" \
        -n -m KR \
        --deepcopy \
        "${MCCOOL}@${RES}" \
        "${OUT}" \
        > "run_balance_${RES}.log" 2>&1
    echo "DONE ${RES}: ${OUT} ($(du -h "${OUT}" | cut -f1))"
done
echo "ALL_BALANCE_DONE"
