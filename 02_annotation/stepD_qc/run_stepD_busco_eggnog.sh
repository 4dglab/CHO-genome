#!/usr/bin/env bash
# PURPOSE: Protein BUSCO + EggNOG functional annotation
# INPUT:   high-confidence proteome
# OUTPUT:  BUSCO / EggNOG results
# TOOL:    BUSCO 6.0.0 / eggNOG-mapper 2.1.12
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# run_stepD_busco_eggnog.sh — Step D: BUSCO dual-set comparison + EggNOG (full genome)
#
# Two-tier annotation strategy (reviewer response):
#   HC  = high-confidence set (Evidence_level L1/L2/L3 + lncRNA), 20,161 genes
#   ALL = high-confidence + L4 ab initio-only predictions,     37,951 genes
# D1a BUSCO protein mode on HC  proteins  -> primary evidence for the manuscript
# D1b BUSCO protein mode on ALL proteins  -> empirical justification that
#         excluding L4 costs little completeness (dual-set comparison)
# D1c BUSCO genome mode on the assembly   -> assembly-level reference context
# D2  EggNOG-mapper functional annotation on HC proteins
set -eo pipefail


BASE=${WORK_DIR}/annotation_final
FINAL=$BASE/stepB4_clean/final
OUT=$BASE/stepD_qc
BUSCO_DL=${BASE%/annotation_final}/busco_downloads
GENOME=${WORK_DIR}/annotation/01_genome/genome.fa
PROT_HC=$FINAL/highconfidence.proteins.fa
PROT_ALL=$FINAL/all.proteins.fa

mkdir -p "$OUT"
cd "$OUT"

echo "[$(date)] D1a: BUSCO protein mode - HIGH-CONFIDENCE set ($(grep -c '^>' $PROT_HC) proteins)"
conda activate busco-env
busco -i "$PROT_HC" -m protein -l mammalia_odb10 -c 32 -o busco_protein_HC \
    --offline --download_path "$BUSCO_DL" 2>&1 | tee D1_protein_HC.log
cp busco_protein_HC/short_summary*.txt busco_protein_HC_summary.txt || true

echo "[$(date)] D1b: BUSCO protein mode - ALL set incl. L4 ($(grep -c '^>' $PROT_ALL) proteins)"
busco -i "$PROT_ALL" -m protein -l mammalia_odb10 -c 32 -o busco_protein_ALL \
    --offline --download_path "$BUSCO_DL" 2>&1 | tee D1_protein_ALL.log
cp busco_protein_ALL/short_summary*.txt busco_protein_ALL_summary.txt || true

echo "[$(date)] D1c: BUSCO genome mode (assembly context, long step)"
if [ ! -s busco_genome/short_summary*.txt ] && [ ! -e busco_genome_summary.txt ]; then
    busco -i "$GENOME" -m genome -l mammalia_odb10 -c 48 -o busco_genome \
        --offline --download_path "$BUSCO_DL" 2>&1 | tee D1_genome.log
    cp busco_genome/short_summary*.txt busco_genome_summary.txt || true
else
    echo "genome mode already done, skipping"
fi

echo "[$(date)] D2: EggNOG-mapper on HIGH-CONFIDENCE proteins"
conda activate annotation-env
emapper.py -i "$PROT_HC" -o eggnog_output --cpu 32 --dbmem --override 2>&1 | tee D2_eggnog.log

echo "[$(date)] ====== BUSCO dual-set comparison summary ======"
for s in busco_protein_HC_summary.txt busco_protein_ALL_summary.txt; do
    echo "--- $s ---"
    grep -E "^# Busco|Complete|Single|Duplicated|Fragmented|Missing|n_markers" "$s" 2>/dev/null | head -8
done
echo ""
head -3 eggnog_output.emapper.annotations 2>/dev/null || true

echo "STEPD_DONE -> $OUT"
