#!/usr/bin/env bash
# PURPOSE: Count onTAD categories (gap / Singleton / level1..level5) at 10 kb.
# INPUT:   ${HI_C_DIR}/TAD/combined/genome.10kb.onTAD.annotation.bed
# OUTPUT:  ${HI_C_DIR}/TAD/combined/onTAD.category_count.10kb.tsv
#          columns: scaffold | category | n_regions | total_bp
#          'Total' row = genome-wide sum over scaffold_1-11
# TOOL:    awk
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

RES=10kb
COMBINED="${HI_C_DIR}/TAD/combined"
ANNOT="${COMBINED}/genome.${RES}.onTAD.annotation.bed"
OUT="${COMBINED}/onTAD.category_count.${RES}.tsv"
BODY="${OUT}.body"

: > "${BODY}"
awk -F'\t' -v OFS='\t' '
  { c[$1 FS $4]++; b[$1 FS $4] += ($3 - $2) }
  END{ for (k in c) { split(k, x, FS); printf "%s\t%s\t%s\t%s\n", x[1], x[2], c[k], b[k] } }
' "${ANNOT}" | sort -V >> "${BODY}"
awk -F'\t' -v OFS='\t' '
  { cg[$4]++; bg[$4] += ($3 - $2) }
  END{
    n=split("gap Singleton level1 level2 level3 level4 level5", L, " ");
    for (i=1; i<=n; i++) { l=L[i]; printf "Total\t%s\t%s\t%s\n", l, cg[l]+0, bg[l]+0 }
  }' "${ANNOT}" >> "${BODY}"

{
  echo "# onTAD category counts, ${RES} resolution (scaffold_1-11)."
  echo "# categories: gap | Singleton | level1..level5 (onTAD nesting depth)."
  echo "# 'Total' row = genome-wide sum over scaffold_1-11."
  printf 'scaffold\tcategory\tn_regions\ttotal_bp\n'
  cat "${BODY}"
} > "${OUT}"
rm -f "${BODY}"
echo "Wrote ${OUT}"
cat "${OUT}"
