#!/usr/bin/env bash
# PURPOSE: Build the genome-wide TAD/gap annotation (Singleton | level1-5 | gap) from
#          the encoded OnTAD table (10 kb).
# INPUT:   ${HI_C_DIR}/TAD/combined/onTAD.10kb.encoded.tsv, ${SCAFFOLD_LENGTHS}
# OUTPUT:  ${HI_C_DIR}/TAD/combined/genome.10kb.onTAD.annotation.bed  (chrom start end label)
# TOOL:    awk / bedtools complement
# env: hic-tools
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   SCAFFOLD_LENGTHS   two-column TSV <scaffold> <length> (cut -f1,2 genome.fa.fai)
#   BEDTOOLS           bedtools executable (hic-tools env)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

RESLABEL=10kb
TAD="${HI_C_DIR}/TAD"
OUTDIR="${TAD}/combined"
LEN="${SCAFFOLD_LENGTHS}"
TSV="${OUTDIR}/onTAD.${RESLABEL}.encoded.tsv"
REGIONS="${OUTDIR}/onTAD.${RESLABEL}.regions.bed"
REGIONS_CLIP="${OUTDIR}/onTAD.${RESLABEL}.regions.clipped.bed"
GAP_RAW="${OUTDIR}/onTAD.${RESLABEL}.gap_raw.bed"
GAP="${OUTDIR}/onTAD.${RESLABEL}.gap.bed"
ANNOT="${OUTDIR}/genome.${RESLABEL}.onTAD.annotation.bed"
mkdir -p "${OUTDIR}"

echo "== [Step 1] Singleton labeling -> 4-col TAD regions =="
tail -n+2 "${TSV}" | awk -F'\t' '
  { if(prev!=""){
      split(prev,p,"\t")
      if(p[1]==$1 && p[6]=="level1" && $6=="level1") printf "%s\t%s\t%s\tSingleton\n", p[1],p[2],p[3]
      else printf "%s\t%s\t%s\t%s\n", p[1],p[2],p[3],p[6]
    }
    if($1!=prev_ch){ prev_ch=$1; prev="" }
    prev=$0 }
  END{ if(prev!=""){ split(prev,p,"\t"); printf "%s\t%s\t%s\t%s\n", p[1],p[2],p[3],p[6] } }
' > "${REGIONS}"

echo "== [Step 2] Clip TAD ends to chromosome length =="
awk -F'\t' 'NR==FNR{len[$1]=$2;next}{ if($3>len[$1])$3=len[$1]; print $1"\t"$2"\t"$3"\t"$4 }' \
  "${LEN}" "${REGIONS}" > "${REGIONS_CLIP}"

echo "== [Step 3] gap = complement of TAD union =="
sort -k1,1 -k2,2n "${REGIONS_CLIP}" > "${OUTDIR}/.sort.regions.${RESLABEL}.tmp"
sort -k1,1 "${LEN}" > "${OUTDIR}/.sort.len.${RESLABEL}.tmp"
"${BEDTOOLS}" complement -i "${OUTDIR}/.sort.regions.${RESLABEL}.tmp" \
  -g "${OUTDIR}/.sort.len.${RESLABEL}.tmp" > "${GAP_RAW}"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"gap"}' "${GAP_RAW}" > "${GAP}"
rm -f "${OUTDIR}/.sort.regions.${RESLABEL}.tmp" "${OUTDIR}/.sort.len.${RESLABEL}.tmp"

echo "== [Step 4] merge gap + TAD regions -> final annotation =="
cat "${GAP}" "${REGIONS_CLIP}" | sort -k1,1 -k2,2n > "${ANNOT}"
echo "annotation lines: $(wc -l < "${ANNOT}")"
echo "DONE -> ${ANNOT}"
