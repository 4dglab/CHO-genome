#!/usr/bin/env bash
# PURPOSE: Step B gene annotation (B0 prep -> B1 Augustus -> B2 EviAnn -> B3 AGAT)
# INPUT:   RNA-seq/Iso-seq + softmasked genome
# OUTPUT:  03_agat/genomic.gff
# TOOL:    HISAT2/StringTie/Augustus/EviAnn/AGAT
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# Master runner for Step B (full genome annotation_final)
# Sequential execution: B0-prep -> B1 -> B2 -> B3 -> B4
# Usage: bash run_stepB_all.sh [--skip-prep]
set -eo pipefail

BASE=${WORK_DIR}/annotation_final
STEP=$BASE/stepB_gene
LOG=$BASE/stepB_gene/logs

mkdir -p "$LOG"

SKIP_PREP=0
[ "${1:-}" = "--skip-prep" ] && SKIP_PREP=1

echo "===== $(date) Step B pipeline start ====="

if [ $SKIP_PREP -eq 0 ]; then
    echo "[$(date)] B0b: RNA-seq HISAT2 + StringTie"
    bash "$STEP/00_prep/run_0b_rnaseq_stringtie.sh" 2>&1 | tee "$LOG/B0b.log"

    echo "[$(date)] B0c: Iso-seq alignment"
    bash "$STEP/00_prep/run_0c_isoseq.sh" 2>&1 | tee "$LOG/B0c.log"

    echo "[$(date)] B0d: build paired.txt"
    bash "$STEP/00_prep/build_paired_txt.sh" 2>&1 | tee "$LOG/B0d.log"
else
    echo "[$(date)] B0-prep skipped (--skip-prep)"
fi

echo "[$(date)] B1: Augustus training + ab initio"
bash "$STEP/01_augustus/run_stepB1_augustus.sh" 2>&1 | tee "$LOG/B1.log"

echo "[$(date)] B2: EviAnn"
bash "$STEP/02_eviann/run_stepB2_eviann.sh" 2>&1 | tee "$LOG/B2.log"

echo "[$(date)] B3: AGAT post-processing"
bash "$STEP/03_agat/run_stepB3_agat.sh" 2>&1 | tee "$LOG/B3.log"

echo "[$(date)] B4: gene filtering/cleaning"
bash "$BASE/stepB4_clean/run_stepB4_clean.sh" 2>&1 | tee "$LOG/B4.log"

echo "===== $(date) Step B pipeline DONE ====="
