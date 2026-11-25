#!/usr/bin/env bash
# PURPOSE: GenomeSyn synteny/repeat visualization
# INPUT:   nucmer results + annotations
# OUTPUT:  GenomeSyn figures
# TOOL:    GenomeSyn 1.2.7
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"


# GenomeSyn rerun with manually-built coords (-cf1).
# Query = scaffold_1-11, gene density -GD1/-GD2 (-GD_win 1Mb), TE tracks,
# NOTE: use conda 'genomesyn' perl; do NOT add localperl PERL5LIB.
set -euo pipefail

conda activate genomesyn
source ${GENOMESYN}/env.sh

GENSYN=${GENOMESYN}/bin/GenomeSyn

WD=${WORK_DIR}/reviewer_response/pur_analysis
RUN="$WD/genomesyn_run_chr"
REF="$RUN/.._ref_chr_contact.fna"
QUERY="$RUN/query_chr1_11.fa"
COORDS="$RUN/ReferencevsQuery1.fixed.delta.filter.coords"
TE1="$RUN/ref_TE.bed"
TE2="$RUN/query_TE.bed"
GD1="$RUN/ref_GD.gff"
GD2="$RUN/query_GD.gff"
PUR2="$RUN/pur_genome1_C1_chr1_11.bed"

cd "$RUN"

echo "=== GenomeSyn -cf1 rerun (scaffold_1-11 + GD + TE + gene2 PUR), start $(date) ==="
"$GENSYN" -t 3 \
  -g1 "$REF" \
  -g2 "$QUERY" \
  -cf1 "$COORDS" \
  -TE1 "$TE1" \
  -TE2 "$TE2" \
  -TE_win 10000 \
  -TE_min 20 \
  -GD1 "$GD1" \
  -GD2 "$GD2" \
  -GD_win 1000000 \
  -gene2 "$PUR2" \
  -c14 "#DC143C" \
  -r Ref \
  -q CHO \
  -pdf 0 2>&1

echo "=== Finished at $(date) ==="
ls -la "$RUN"/*.svg 2>/dev/null || ls -la "$RUN/"
