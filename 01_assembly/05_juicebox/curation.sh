#!/usr/bin/env bash
# PURPOSE: juicer post after manual curation (review4)
# INPUT:   cho_yahs_v3.JBAT.review4.assembly + liftover.agp + purged_clean.fa
# OUTPUT:  cho_v3_review4.fa + .agp
# TOOL:    juicer post
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -e
conda activate yahs-env

WORKDIR="${WORK_DIR}/yahs/output_v3"
ASM_FA="${WORK_DIR}/purge_dups_v3/purged_clean.fa"
LIFTOVER_AGP="${WORKDIR}/cho_yahs_v3.JBAT.liftover.agp"
REVIEWED="${WORKDIR}/cho_yahs_v3.JBAT.review4.assembly"
OUTDIR="${WORKDIR}/juicer_post_review4"
FINAL_PREFIX="${OUTDIR}/cho_v3_review4"

mkdir -p "${OUTDIR}"

echo "===== $(date) juicer post v3 review4 ====="
echo "Input contigs:   ${ASM_FA}"
echo "Liftover AGP:    ${LIFTOVER_AGP}"
echo "Reviewed RSM:    ${REVIEWED}"
echo "Output prefix:   ${FINAL_PREFIX}"

juicer post -o "${FINAL_PREFIX}" "${REVIEWED}" "${LIFTOVER_AGP}" "${ASM_FA}"

if [ -s "${FINAL_PREFIX}.FINAL.fa" ]; then
    mv "${FINAL_PREFIX}.FINAL.fa" "${FINAL_PREFIX}.fa"
    mv "${FINAL_PREFIX}.FINAL.agp" "${FINAL_PREFIX}.agp"
fi

echo "===== $(date) Done ====="
awk '/^>/{if(n++) print len; len=0; next} {len+=length($0)} END{print len}' "${FINAL_PREFIX}.fa" | \
    sort -rn | awk '{vals[NR]=$1; sum+=$1}
    END{
        half=sum/2; cum=0;
        for(i=1;i<=NR;i++){cum+=vals[i]; if(cum>=half && !found){printf "N50: %.2f Mb\n", vals[i]/1e6; found=1}}
        printf "Longest: %.2f Mb  Total: %.2f Gb  Scaffolds: %d\n", vals[1]/1e6, sum/1e9, NR
    }'
echo -n "N count: "; tr -dc 'Nn' < "${FINAL_PREFIX}.fa" | wc -c
