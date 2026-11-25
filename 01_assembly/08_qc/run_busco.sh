#!/usr/bin/env bash
# PURPOSE: BUSCO genome completeness assessment
# INPUT:   cho_review4_gapfilled_polished.fa
# OUTPUT:  busco_polished/ results
# TOOL:    BUSCO 6.0.0 (mammalia_odb10)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -eo pipefail
conda activate busco-env

POLISHED="${WORK_DIR}/nextpolish2_review4_gapfilled/output/cho_review4_gapfilled_polished.fa"
OUTDIR="${WORK_DIR}/qc_final/busco"
THREADS=32

echo "========== $(date) Step 7 QC: BUSCO on polished genome =========="
echo "Genome:  ${POLISHED}"
echo "Threads: ${THREADS}"
echo "Lineage: mammalia_odb10"

busco -i "${POLISHED}" \
    -l mammalia_odb10 \
    -m genome \
    -c ${THREADS} \
    -o busco_polished \
    --out_path "${OUTDIR}" \
    --download_path ${WORK_DIR}/busco/busco_downloads \
    --offline \
    --force \
    2>&1 | tee "${OUTDIR}/busco.log"

echo ""
echo "========== $(date) BUSCO DONE =========="
head -20 "${OUTDIR}/busco_polished/short_summary.txt"
