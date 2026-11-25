#!/usr/bin/env bash
# PURPOSE: Merge HiTE+Dfam repeat annotation, stats, softmask
# INPUT:   HiTE + Dfam GFF3
# OUTPUT:  merge_all.repeat.gff3 + genome.softmasked.fa
# TOOL:    custom scripts
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# Step A2-A5: merge + dedup -> stats -> softmask (run AFTER HiTE finishes)
set -eo pipefail

BASE=${WORK_DIR}/annotation_final
SCRIPTS=${BASE}/scripts
REPEAT_DIR=${BASE}/stepA_repeat
HITE_OUT=${BASE}/hite_out/HiTE.out
DFAM_OUT=${WORK_DIR}/annotation/04_repeatmasker/genome.fa.out
GENOME=${WORK_DIR}/annotation/01_genome/genome.fa
GENOME_SIZE=2670329556

mkdir -p "$REPEAT_DIR"

echo "[$(date)] A3: merge HiTE + Dfam -> merge_all.repeat.gff3"
python3 "$SCRIPTS/merge_repeat_gff3.py" \
    --hite "$HITE_OUT" \
    --dfam "$DFAM_OUT" \
    --genome "$GENOME" \
    --out "$REPEAT_DIR/merge_all.repeat.gff3" 2>&1 | tee "$REPEAT_DIR/merge_stats.log"

echo "[$(date)] A4: repeat content stats"
python3 "$SCRIPTS/repeat_content_stats.py" \
    --gff "$REPEAT_DIR/merge_all.repeat.gff3" \
    --genome "$GENOME_SIZE" \
    --outdir "$REPEAT_DIR" 2>&1 | tee "$REPEAT_DIR/repeat_stats.log"

echo "[$(date)] A5: soft-mask genome for Augustus"
python3 "$SCRIPTS/softmask_genome.py" \
    "$GENOME" \
    "$REPEAT_DIR/merge_all.repeat.gff3" \
    "$BASE/genome.softmasked.fa" 2>&1 | tee "$REPEAT_DIR/softmask.log"

echo "[$(date)] Step A (whole-genome repeat) done"
