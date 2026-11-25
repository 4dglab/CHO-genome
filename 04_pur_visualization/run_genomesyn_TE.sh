#!/usr/bin/env bash
# PURPOSE: GenomeSyn synteny/repeat visualization
# INPUT:   nucmer results + annotations
# OUTPUT:  GenomeSyn figures
# TOOL:    GenomeSyn 1.2.7
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"


# Run classic GenomeSyn 1.x for repeat-content (TE) collinearity plot
# Reference = 11-chr CriGri (chr_contact_chr1), Query = CHO assembly
# GenomeSyn calls MUMmer automatically (no -cf1)
set -euo pipefail

conda activate genomesyn
source ${GENOMESYN}/env.sh

GENSYN=${GENOMESYN}/bin/GenomeSyn

REF="${RAW_DIR}/genome_assembly/assembly_ref/GCF_003668045.3_CriGri-PICRH-1.0_genomic_chr_contact_chr1.fna"
QUERY="${WORK_DIR}/nextpolish2_review4_gapfilled/block46_translocation/output/cho_review5_gapfilled_polished.fa"
TE_REF="${WORK_DIR}/reviewer_response/pur_analysis/genomesyn_te/ref_TE_genomesyn.bed"
TE_QUERY="${WORK_DIR}/reviewer_response/pur_analysis/genomesyn_te/query_TE_genomesyn.bed"

RUN="${WORK_DIR}/reviewer_response/pur_analysis/genomesyn_run"
mkdir -p "$RUN"
cd "$RUN"

echo "=== Running GenomeSyn (repeat-content synteny, MUMmer auto-alignment) ==="
"$GENSYN" -t 3 \
  -g1 "$REF" \
  -g2 "$QUERY" \
  -TE1 "$TE_REF" \
  -TE2 "$TE_QUERY" \
  -TE_win 10000 \
  -TE_min 20 \
  -r Ref \
  -q CHO \
  -pdf 0 2>&1

echo "=== Finished, output files: ==="
ls -la "$RUN"/*.svg 2>/dev/null || ls -la "$RUN/"
