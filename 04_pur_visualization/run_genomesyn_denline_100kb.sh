#!/usr/bin/env bash
# PURPOSE: GenomeSyn synteny/repeat visualization
# INPUT:   nucmer results + annotations
# OUTPUT:  GenomeSyn figures
# TOOL:    GenomeSyn 1.2.7
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"


# GenomeSyn rerun with:
#  - TE as DENSITY LINES (not per-window rects): via patched copy GenomeSyn_densityline
#    (source line 1138 $snp_TE forced to 0 -> black TE density polylines)
#  - Gene density color -> #A6CEE3 (-c13)
#  - gene2 PUR track on CHO (query) side (-gene2, -c14 #DC143C)
# Query = scaffold_1-11, -cf1 manually-built coords (space-in-reffpath fix), -GD_win 1Mb.
# NOTE: use conda 'genomesyn' perl; do NOT add localperl PERL5LIB.
set -euo pipefail

conda activate genomesyn
source ${GENOMESYN}/env.sh

GENSYN=${GENOMESYN}/bin/GenomeSyn_densityline

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

echo "=== GenomeSyn_densityline -cf1 (TE density line + GD#A6CEE3 + gene2 PUR), start $(date) ==="
"$GENSYN" -t 3 \
  -g1 "$REF" \
  -g2 "$QUERY" \
  -cf1 "$COORDS" \
  -TE1 "$TE1" \
  -TE2 "$TE2" \
  -TE_win 100000 \
  -GD1 "$GD1" \
  -GD2 "$GD2" \
  -GD_win 1000000 \
  -c13 "#A6CEE3" \
  -gene2 "$PUR2" \
  -c14 "#DC143C" \
  -r Ref \
  -q CHO \
  -pdf 0 2>&1

echo "=== Finished at $(date) ==="
ls -la "$RUN"/*.svg 2>/dev/null || ls -la "$RUN/"
