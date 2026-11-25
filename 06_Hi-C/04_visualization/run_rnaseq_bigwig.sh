#!/usr/bin/env bash
# PURPOSE: Align A485 RNA-seq (control Cro and treated 10uM) and generate RPKM bigwigs
#          used as the RNA-seq tracks in the Hi-C/epigenetics figures.
# INPUT:   ${A485_RNA_DIR}/{Cro-1,Cro-2,Cro-3,10uM-1,10uM-2,10uM-3}/*_R{1,2}.fq.gz
# OUTPUT:  ${HI_C_DIR}/bigwig/Cro.sorted.newref.bigwig
#          ${HI_C_DIR}/bigwig/10uM.sorted.newref.bigwig
# TOOL:    HISAT2 2.2.3, samtools, deepTools bamCoverage 3.5.6 (--binSize 10 --normalizeUsing RPKM)
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   GENOME_FA          reference genome fasta with .fai index
#   GENOME_GFF         gene annotation GFF3
#   A485_RNA_DIR       clean A485 reads <sample>/<sample>_R{1,2}.fq.gz
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

OUT="${HI_C_DIR}/bigwig"
IDX="${HI_C_DIR}/index/gapfill_v4_ht2"
SPLICE="${OUT}/gapfill_v4.splice"
CLEAN="${A485_RNA_DIR}"
T="${THREADS}"
mkdir -p "${OUT}" "${HI_C_DIR}/index"

declare -A GROUPS=(
  [Cro]="Cro-1 Cro-2 Cro-3"
  [10uM]="10uM-1 10uM-2 10uM-3"
)

# ---- align with HISAT2 + samtools (annotation-env) ----
eval "$(conda shell.bash hook)" 2>/dev/null || true
conda activate annotation-env

if [ ! -s "${IDX}.1.ht2" ]; then
    echo "[$(date)] hisat2-build"
    hisat2-build --threads "${T}" "${GENOME_FA}" "${IDX}"
fi
if [ -s "${GENOME_GFF}" ]; then
    hisat2_extract_splice_sites.py "${GENOME_GFF}" > "${SPLICE}" 2>/dev/null || true
fi
SPLICE_ARG=()
[ -s "${SPLICE}" ] && SPLICE_ARG=(--known-splicesite-infile "${SPLICE}")

for G in Cro 10uM; do
    for s in ${GROUPS[$G]}; do
        echo "[$(date)] hisat2 ${s}"
        hisat2 --phred33 -x "${IDX}" \
            -1 "${CLEAN}/${s}/${s}_R1.fq.gz" -2 "${CLEAN}/${s}/${s}_R2.fq.gz" \
            -p "${T}" "${SPLICE_ARG[@]}" 2> "${OUT}/${s}.ht2.log" \
          | samtools view -h -b -f 3 -F 12 -F 256 -q 20 -o "${OUT}/${s}.unsorted.bam" -
        samtools sort -@ "${T}" -o "${OUT}/${s}.sort.bam" "${OUT}/${s}.unsorted.bam"
        rm -f "${OUT}/${s}.unsorted.bam"
        samtools index -c "${OUT}/${s}.sort.bam"
    done
    echo "[$(date)] merge ${G}"
    samtools merge -@ "${T}" -o "${OUT}/${G}.merged.sorted.newref.bam" \
        $(for s in ${GROUPS[$G]}; do echo "${OUT}/${s}.sort.bam"; done)
    samtools index -c "${OUT}/${G}.merged.sorted.newref.bam"
done

# ---- bamCoverage (deeptools-env) ----
conda activate deeptools-env
for G in Cro 10uM; do
    echo "[$(date)] bamCoverage ${G} (RPKM, bin10)"
    bamCoverage -b "${OUT}/${G}.merged.sorted.newref.bam" \
        -o "${OUT}/${G}.sorted.newref.bigwig" \
        --binSize 10 --normalizeUsing RPKM -p "${T}"
done
echo "[$(date)] ALL_RNA_BIGWIG_DONE"
ls -lh "${OUT}"/*.sorted.newref.bigwig
