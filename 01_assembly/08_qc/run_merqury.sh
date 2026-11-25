#!/usr/bin/env bash
# PURPOSE: Merqury QV / k-mer completeness assessment
# INPUT:   cho_review4_gapfilled_polished.fa + NGS clean reads
# OUTPUT:  polished_qv.* results
# TOOL:    Merqury 1.4.1 / meryl
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -eo pipefail

POLISHED="${WORK_DIR}/nextpolish2_review4_gapfilled/output/cho_review4_gapfilled_polished.fa"
NGS_R1="${WORK_DIR}/nextpolish2/ngs_clean/R1.clean.fq.gz"
NGS_R2="${WORK_DIR}/nextpolish2/ngs_clean/R2.clean.fq.gz"
OUTDIR="${WORK_DIR}/qc_final/merqury"
THREADS=48

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

LOG="${OUTDIR}/run.log"
exec > >(tee -a "${LOG}") 2>&1

echo "========== $(date) Merqury QV analysis =========="
echo "Genome:   ${POLISHED}"
echo "NGS R1:   ${NGS_R1} ($(du -h ${NGS_R1} | cut -f1))"
echo "NGS R2:   ${NGS_R2} ($(du -h ${NGS_R2} | cut -f1))"
echo "Threads:  ${THREADS}"

# ============================================================
# Step 1: Install merqury
# ============================================================
echo "[$(date)] Step 1: Installing merqury..."
conda activate nextpolish2-env
conda install -y -n nextpolish2-env merqury=1.4.1 -c bioconda -c conda-forge 2>&1 | tail -5
echo "[$(date)] merqury installed"

# ============================================================
# Step 2: Build genome meryl k=21 DB
# ============================================================
echo "[$(date)] Step 2: Building genome meryl k=21 DB..."
GENOME_MERYL="${OUTDIR}/genome.k21.meryl"
meryl count k=21 threads=${THREADS} output "${GENOME_MERYL}" "${POLISHED}" 2>&1 | tail -5
echo "[$(date)] Genome meryl DONE"
du -sh "${GENOME_MERYL}"

# ============================================================
# Step 3: Build NGS reads meryl k=21 DB
# ============================================================
echo "[$(date)] Step 3: Building NGS meryl k=21 DB..."
NGS_MERYL="${OUTDIR}/ngs.k21.meryl"
meryl count k=21 threads=${THREADS} memory=200g \
    output "${NGS_MERYL}" "${NGS_R1}" "${NGS_R2}" 2>&1 | tail -5
echo "[$(date)] NGS meryl DONE"
du -sh "${NGS_MERYL}"

# ============================================================
# Step 4: Run merqury
# ============================================================
echo "[$(date)] Step 4: Running merqury.sh..."
export MERQURY=${HOME}/anaconda3/envs/nextpolish2-env/share/merqury
Rscript ${HOME}/anaconda3/envs/nextpolish2-env/lib/R/library/argparse/R/parse.R 2>/dev/null || true

merqury.sh "${NGS_MERYL}" "${GENOME_MERYL}" polished_qv 2>&1 | tail -20

echo ""
echo "========== $(date) Merqury DONE =========="
echo "Output directory: ${OUTDIR}"
ls -lh polished_qv.* 2>/dev/null

# Show key metrics
if [ -f polished_qv.qv ]; then
    echo ""
    echo "=== QV Report ==="
    cat polished_qv.qv
fi
