#!/usr/bin/env bash
# PURPOSE: ncRNA annotation (Rfam cmscan + tRNAscan-SE)
# INPUT:   genome.fa
# OUTPUT:  ncRNA.gff / tRNA.gff
# TOOL:    Infernal cmscan / tRNAscan-SE 2.0
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# run_stepC_ncrna.sh — Step C: ncRNA annotation (FULL genome, 2.67 Gb / 74 contigs)
#   C1. tRNAscan-SE v2.0.13 -> tRNA.gff          (--thread 72)
#   C2. Infernal cmscan vs Rfam -> ncRNA.gff     (--cpu 72)
#   C3. extract rRNA/snRNA/miRNA sub-files (GFF3 header preserved)
#
# Input genome: UNMASKED (01_genome/genome.fa) per methodology.
# Serial execution; idempotent resume guards on each step.
set -eo pipefail

conda activate annotation-env

BASE=${WORK_DIR}/annotation_final
GENOME=${WORK_DIR}/annotation/01_genome/genome.fa
RFAM=${WORK_DIR}/databases/rfam/Rfam.cm
OUT=$BASE/stepC_ncrna
SCRIPTS=$BASE/scripts
THREADS=72

mkdir -p "$OUT"
cd "$OUT"

echo "========== $(date) Step C: ncRNA annotation (genome=unmasked, threads=$THREADS) =========="

# ---------- C1: tRNAscan-SE ----------
if [ ! -s tRNA.gff ]; then
    echo "[$(date)] C1: tRNAscan-SE v2.0.13 (-E -X 20 -z 8 --thread $THREADS)"
    tRNAscan-SE -E -X 20 -z 8 --thread $THREADS \
        -j tRNA.gff -o tRNA.out -m tRNA.stats "$GENOME" 2>&1 | tee C1_tRNAscan.log
    head -2 tRNA.gff
else
    echo "[$(date)] C1: tRNA.gff exists, skipping"
fi

# ---------- C2: cmscan vs Rfam ----------
if [ ! -s ncRNA.gff ]; then
    echo "[$(date)] C2: cmscan vs Rfam (E<0.001, --cpu $THREADS)"
    cmscan --rfam --fmt 2 -E 0.001 \
        --tblout rfam.tblout --cpu $THREADS \
        "$RFAM" "$GENOME" > cmscan.out 2> cmscan.err
    python3 "$SCRIPTS/cmscan_tblout2gff3.py" \
        --tblout rfam.tblout --genome "$GENOME" \
        --evalue 0.001 --out ncRNA.gff
else
    echo "[$(date)] C2: ncRNA.gff exists, skipping"
fi

# ---------- C3: extract sub-files (GFF3 header preserved) ----------
echo "[$(date)] C3: extract rRNA/snRNA/miRNA sub-files"
echo '##gff-version 3' > rRNA.gff
grep -P '\trRNA\t' ncRNA.gff >> rRNA.gff || true
echo '##gff-version 3' > snRNA.gff
grep -P '\tsnRNA\t' ncRNA.gff >> snRNA.gff || true
echo '##gff-version 3' > miRNA.gff
grep -P '\tmiRNA\t' ncRNA.gff >> miRNA.gff || true

echo "[$(date)] summary:"
for f in ncRNA.gff rRNA.gff snRNA.gff miRNA.gff tRNA.gff; do
    printf "%-12s features=%s\n" "$f" "$(grep -vc '^#' "$f")"
done

echo "STEPC_DONE -> $OUT"
