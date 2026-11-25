#!/usr/bin/env bash
# PURPOSE: GenomeSyn synteny/repeat visualization
# INPUT:   nucmer results + annotations
# OUTPUT:  GenomeSyn figures
# TOOL:    GenomeSyn 1.2.7
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"


# Rerun GenomeSyn with:
#  - query = extracted scaffold_1-11 (name corresponds to ref chr1-10 + X)
#  - gene density tracks: -GD1 (ref GFF) -GD2 (query GFF) -GD_win 1000000
#  - TE tracks: -TE1 -TE2 (cleaned BED, header stripped)
#  - note: use conda 'genomesyn' env perl which already has SVG.pm (do NOT
#    add localperl PERL5LIB -- it is build on another host and version-mismatched)
# GenomeSyn calls MUMmer automatically (no -cf1)
set -euo pipefail

conda activate genomesyn
source ${GENOMESYN}/env.sh

GENSYN=${GENOMESYN}/bin/GenomeSyn

WD=${WORK_DIR}/reviewer_response/pur_analysis
RUN="$WD/genomesyn_run_chr"
REF="${RAW_DIR}/genome_assembly/assembly_ref/GCF_003668045.3_CriGri-PICRH-1.0_genomic_chr_contact_chr1.fna"
QUERY="$RUN/query_chr1_11.fa"
TE1="$RUN/ref_TE.bed"
TE2="$RUN/query_TE.bed"
GD1="$RUN/ref_GD.gff"
GD2="$RUN/query_GD.gff"

cd "$RUN"

echo "=== GenomeSyn rerun (scaffold_1-11 + GD + TE), start $(date) ==="
"$GENSYN" -t 3 \
  -g1 "$REF" \
  -g2 "$QUERY" \
  -TE1 "$TE1" \
  -TE2 "$TE2" \
  -TE_win 10000 \
  -TE_min 20 \
  -GD1 "$GD1" \
  -GD2 "$GD2" \
  -GD_win 1000000 \
  -r Ref \
  -q CHO \
  -pdf 0 2>&1

echo "=== Finished at $(date) ==="
ls -la "$RUN"/*.svg 2>/dev/null || ls -la "$RUN/"
