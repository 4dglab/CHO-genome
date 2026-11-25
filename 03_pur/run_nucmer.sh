#!/usr/bin/env bash
# PURPOSE: nucmer whole-genome alignment (PUR / synteny)
# INPUT:   final genome + CriGri-PICRH-1.0 reference
# OUTPUT:  ref_vs_cho.delta + coords
# TOOL:    MUMmer4 nucmer
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"


set -eo pipefail

conda activate mummer-env

export PATH="${GNUPLOT_BIN}:$PATH"
export GNUPLOT_PS_DIR="${GNUPLOT_PS_DIR}"

WORKDIR="${WORK_DIR}/collinearity_v3_review5"
REF="${RAW_DIR}/genome_assembly/assembly ref/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
QRY="${WORK_DIR}/yahs/output_v3/juicer_post_review5/cho_v3_review5.fa"
PREFIX="GCF_CriGri_PICRH_vs_cho_v3_review5"
THREADS=48

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

echo ">>>> $(date) nucmer v3 review5 start <<<<"
echo "Ref:   ${REF}"
echo "Query: ${QRY}"

echo "[$(date)] Step 1: nucmer..."
nucmer --maxgap=500 --mincluster=100 -t ${THREADS} -p "${PREFIX}" "${REF}" "${QRY}" 2> logs/nucmer.log
echo "[$(date)] done: $(du -h ${PREFIX}.delta | cut -f1)"

echo "[$(date)] Step 2: delta-filter..."
delta-filter -m -l 2000 "${PREFIX}.delta" > "${PREFIX}.filtered.delta" 2> logs/delta_filter.log
echo "[$(date)] done: $(du -h ${PREFIX}.filtered.delta | cut -f1)"

echo "[$(date)] Step 3: show-coords..."
show-coords -c -l -T "${PREFIX}.filtered.delta" > "${PREFIX}.coords.txt" 2> logs/show_coords.log
echo "[$(date)] done: $(tail -n +5 ${PREFIX}.coords.txt | wc -l) alignments"

echo "[$(date)] Step 4: Chr assignment..."
awk -F'\t' '
BEGIN { print "scaffold\tchromosome\ttotal_aligned_bp" }
NR<=4 {next}
{ref=$12; qry=$13; alen=$5+0; bases[qry SUBSEP ref]+=alen; scaffolds[qry]=1; chrs[ref]=1}
END {
    for(s in scaffolds){max=0; best="unassigned"; total=0
        for(r in chrs){b=bases[s SUBSEP r]+0; if(b>0){total+=b; if(b>max){max=b; best=r}}}
        printf "%s\t%s\t%d\n", s, best, total}
}' "${PREFIX}.coords.txt" > "${PREFIX}.chr_assignment.tsv"

echo "[$(date)] Assigned: $(awk -F'\t' '$2!="unassigned"' ${PREFIX}.chr_assignment.tsv | tail -n +2 | wc -l), Unassigned: $(awk -F'\t' '$2=="unassigned"' ${PREFIX}.chr_assignment.tsv | wc -l)"

echo "[$(date)] Step 5: mummerplot..."
mummerplot -postscript -p "${PREFIX}" "${PREFIX}.filtered.delta" 2> logs/mummerplot.log
echo "[$(date)] done: ${PREFIX}.ps"

if [ -f "${PREFIX}.ps" ]; then
    conda activate && python3 -c "
from PIL import Image
Image.open('${WORKDIR}/${PREFIX}.ps').save('${WORKDIR}/${PREFIX}.png')
" 2>/dev/null && echo "[$(date)] png saved"
fi

echo ">>>> $(date) COMPLETE <<<<"
ls -lh "${PREFIX}.delta" "${PREFIX}.filtered.delta" "${PREFIX}.coords.txt" "${PREFIX}.chr_assignment.tsv" "${PREFIX}.ps" "${PREFIX}.png" 2>/dev/null

echo ""
echo "=== Per-chromosome ==="
awk -F'\t' 'NR>1{chr=$2; if(chr=="unassigned")chr="UNASSIGNED"; count[chr]++; bases[chr]+=$3}
END{for(c in count) printf "%-15s  %2d scaffolds  %12d bp\n", c, count[c], bases[c]}' "${PREFIX}.chr_assignment.tsv" | sort
