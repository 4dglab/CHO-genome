#!/usr/bin/env bash
# PURPOSE: HiTE de novo transposon prediction
# INPUT:   genome fasta
# OUTPUT:  HiTE output (TE library/annotation)
# TOOL:    HiTE 3.3.3
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# Step A1: HiTE de novo repeat annotation on the FULL genome (2.67 Gb, 74 contigs)
set -eo pipefail

conda activate HiTE

HITE=${WORK_DIR}/HiTE-master
BASE=${WORK_DIR}/annotation_final
GENOME=${WORK_DIR}/annotation/01_genome/genome.fa
OUT_DIR=${BASE}/hite_out
WORK_DIR=${BASE}/hite_work

export PATH="$HITE/module:$HITE/tools:$PATH"

mkdir -p "$OUT_DIR" "$WORK_DIR"

# do NOT cd into OUT_DIR (HiTE clears it and getcwd() fails)
cd "$BASE"
python "$HITE/main.py" \
    --genome "$GENOME" \
    --thread 64 \
    --annotate 1 \
    --out_dir "$OUT_DIR" \
    --work_dir "$WORK_DIR" \
    2>&1 | tee "$OUT_DIR/hite_run.log"
