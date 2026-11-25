#!/usr/bin/env bash
# PURPOSE: Replace the color column of OnTAD .bed with levelN labels (level1..level5).
# INPUT:   ${HI_C_DIR}/TAD/10kb.onTAD/onTAD.K1.scaffold_*.bed
# OUTPUT:  same .bed files with col9 = levelN (originals backed up to *.bed.backup)
# TOOL:    awk
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

RESLABEL=10kb
OUTDIR="${HI_C_DIR}/TAD/${RESLABEL}.onTAD"

for bed in "${OUTDIR}"/onTAD.K1.scaffold_*.bed; do
    [ -e "$bed" ] || continue
    [ -e "${bed}.backup" ] || cp "$bed" "${bed}.backup"
    awk -v OFS='\t' '
      /^track/{print; next}
      {
        c=$9
        if      (c=="56,108,176"){$9="level1"}
        else if (c=="127,201,127"){$9="level2"}
        else if (c=="190,174,212"){$9="level3"}
        else if (c=="253,192,134"){$9="level4"}
        else if (c=="255,0,0"){$9="level5"}
        else {$9="level"c}
        print
      }' "$bed" > "${bed}.tmp"
    mv "${bed}.tmp" "$bed"
    echo "LABELED: $bed"
done
echo "ALL_LABELED"
