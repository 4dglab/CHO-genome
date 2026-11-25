#!/usr/bin/env bash
# PURPOSE: tidk telomere repeat (TTAGGG) identification
# INPUT:   cho_review4_gapfilled_polished.fa
# OUTPUT:  cho_telo_telomeric_repeat_windows.tsv
# TOOL:    tidk 0.2.65
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -eo pipefail
conda activate hifiasm-env

POLISHED="${WORK_DIR}/nextpolish2_review4_gapfilled/output/cho_review4_gapfilled_polished.fa"
OUTDIR="${WORK_DIR}/qc_final/tidk"
THREADS=16

echo "========== $(date) Step 7 QC: tidk telomere identification =========="
echo "Genome:  ${POLISHED}"
echo "Threads: ${THREADS}"
echo "Repeat:  TTAGGG (vertebrate telomere)"

mkdir -p "${OUTDIR}"

cd "${OUTDIR}"
tidk search \
    --string TTAGGG \
    --output cho_telo \
    --dir . \
    "${POLISHED}" \
    2>&1 | tee "${OUTDIR}/tidk_search.log"

echo ""
echo "========== $(date) tidk DONE =========="

# Summarize
echo ""
echo "=== Telomere arm statistics ==="
NLINES=$(wc -l < "${OUTDIR}/cho_telo_telomeric_repeat_windows.tsv")
echo "Total telomeric repeat windows: $((NLINES - 1))"

# Count scaffolds with telomeric repeats at either end
python3 << 'PYEOF'
import sys
tsv = "${WORK_DIR}/qc_final/tidk/cho_telo_telomeric_repeat_windows.tsv"
ends_with_telomere = set()
with open(tsv) as f:
    header = f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) < 3:
            continue
        scaffold = cols[0]
        start = int(cols[1])
        # telomere at start (within first 500 bp) or end
        if start < 500:
            ends_with_telomere.add(f"{scaffold}_5p")
        # can't determine end position from TSV alone; approximate from total windows per scaffold
        # for now just count 5' telomeres
n_scaffolds = 74
print(f"Scaffolds with 5'-end telomere (TTAGGG within 500bp): {len(ends_with_telomere)}")
PYEOF
