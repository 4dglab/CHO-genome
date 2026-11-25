#!/usr/bin/env bash
# PURPOSE: Identify hierarchical (nested) TADs with OnTAD at 10 kb, per scaffold.
# INPUT:   ${HIC_FILE}, ${SCAFFOLD_LENGTHS}
# OUTPUT:  ${HI_C_DIR}/TAD/10kb.onTAD/onTAD.K1.<scaffold>.{tad,bed}
# TOOL:    OnTAD v1.4 (-bedout)
# env: hic-tools
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   HIC_FILE           merged .hic matrix (HiC-Pro valid pairs -> juicer pre)
#   SCAFFOLD_LENGTHS   two-column TSV <scaffold> <length> (cut -f1,2 genome.fa.fai)
#   ONTAD_BIN          OnTAD v1.4 executable (needs java, hic-tools env)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

RES=10000
RESLABEL=10kb
OUTDIR="${HI_C_DIR}/TAD/${RESLABEL}.onTAD"
mkdir -p "${OUTDIR}"

echo "===== $(date) OnTAD ${RESLABEL} ====="
while read -r scaffold len; do
    [ -n "${scaffold}" ] || continue
    prefix="${OUTDIR}/onTAD.K1.${scaffold}"
    if [ -s "${prefix}.tad" ] && [ -s "${prefix}.bed" ]; then
        echo "SKIP (exists): ${prefix}"
        continue
    fi
    echo "RUN ${scaffold} (len=${len})"
    "${ONTAD_BIN}" "${HIC_FILE}" -bedout "${scaffold}" "${len}" "${RES}" \
        -o "${prefix}" > "${OUTDIR}/log.${scaffold}.txt" 2>&1
done < "${SCAFFOLD_LENGTHS}"

echo "===== $(date) ALL_ONTAD_DONE ====="
