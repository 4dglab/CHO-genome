#!/usr/bin/env bash
# PURPOSE: YaHS Hi-C chromosome-scale scaffolding
# INPUT:   purged_clean.fa + Hi-C clean reads
# OUTPUT:  cho_yahs_v3_scaffolds_final.fa + .agp
# TOOL:    YaHS 1.2.2 / BWA-MEM
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -e

conda activate yahs-env

# ============================================================
# YaHS Hi-C Scaffolding — on cleaned purged contigs
# Input: purged_clean.fa (119 contigs, 2.67 Gb, 0 non-ACGT)
# ============================================================

WORKDIR="${WORK_DIR}/yahs"
OUTDIR="${WORKDIR}/output_v3"
ASSEMBLY="${WORK_DIR}/purge_dups_v3/purged_clean.fa"

HIC_R1_1="${WORK_DIR}/hi-c_data/clean/CH6F_G3-1-1_R1.fq.gz"
HIC_R1_2="${WORK_DIR}/hi-c_data/clean/CH6F_G3-1-1_R2.fq.gz"
HIC_R2_1="${WORK_DIR}/hi-c_data/clean/CH6F_G3-1-2_R1.fq.gz"
HIC_R2_2="${WORK_DIR}/hi-c_data/clean/CH6F_G3-1-2_R2.fq.gz"

THREADS=48
PREFIX="cho_yahs_v3"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

LOG="${OUTDIR}/run.log"
exec > >(tee -a "${LOG}") 2>&1

echo "===== $(date) YaHS v3 Hi-C Scaffolding (cleaned) ====="
echo "Assembly: ${ASSEMBLY}"
echo "Prefix:   ${PREFIX}"
echo "Enzyme:   GATC (DpnII)"
echo "Threads:  ${THREADS}"
echo ""

N_CTG=$(grep -c "^>" "${ASSEMBLY}")
TOTAL_BP=$(grep -v "^>" "${ASSEMBLY}" | tr -d '\n' | wc -c)
echo "[$(date)] Input: ${N_CTG} contigs, ${TOTAL_BP} bp"

for f in "${ASSEMBLY}" "${HIC_R1_1}" "${HIC_R1_2}" "${HIC_R2_1}" "${HIC_R2_2}"; do
    if [ ! -f "$f" ]; then
        echo "[$(date)] ERROR: Missing input: $f"
        exit 1
    fi
done
echo "[$(date)] All inputs verified"

# ============================================================
# Step 1: BWA index
# ============================================================
echo "===== $(date) Step 1: BWA index ====="

BWA_PREFIX="$(dirname "${ASSEMBLY}")/$(basename "${ASSEMBLY}")"

if [ -f "${BWA_PREFIX}.bwt" ] && [ -f "${BWA_PREFIX}.amb" ]; then
    echo "[$(date)] BWA index exists, skip"
else
    cd "$(dirname "${ASSEMBLY}")"
    bwa index -p "$(basename "${ASSEMBLY}")" "$(basename "${ASSEMBLY}")"
fi

echo "===== $(date) Step 1 done ====="

# ============================================================
# Step 2: BWA mem Hi-C → contigs
# ============================================================
echo "===== $(date) Step 2a: BWA mem sample 1 ====="

bwa mem -5SP -t ${THREADS} \
    "${ASSEMBLY}" \
    "${HIC_R1_1}" "${HIC_R1_2}" \
    | samtools view -@ 16 -bS - \
    | samtools sort -@ 16 -o "${OUTDIR}/hic_sample1.bam" -

echo "[$(date)] Sample 1 done"

echo "===== $(date) Step 2b: BWA mem sample 2 ====="

bwa mem -5SP -t ${THREADS} \
    "${ASSEMBLY}" \
    "${HIC_R2_1}" "${HIC_R2_2}" \
    | samtools view -@ 16 -bS - \
    | samtools sort -@ 16 -o "${OUTDIR}/hic_sample2.bam" -

echo "[$(date)] Sample 2 done"

echo "===== $(date) Step 2c: Merge BAMs ====="

samtools merge -@ 16 "${OUTDIR}/hic_aligned.bam" \
    "${OUTDIR}/hic_sample1.bam" "${OUTDIR}/hic_sample2.bam"
samtools index -@ 16 "${OUTDIR}/hic_aligned.bam"

rm "${OUTDIR}/hic_sample1.bam" "${OUTDIR}/hic_sample2.bam"

BAM_SIZE=$(du -h "${OUTDIR}/hic_aligned.bam" | cut -f1)
echo "[$(date)] Merged BAM: ${BAM_SIZE}"

# ============================================================
# Step 3: YaHS scaffolding
# ============================================================
echo "===== $(date) Step 3: YaHS scaffolding (DpnII GATC) ====="

yahs \
    -e GATC \
    -l 10000 \
    -o "${PREFIX}" \
    "${ASSEMBLY}" \
    "${OUTDIR}/hic_aligned.bam"

echo "===== $(date) Step 3 done ====="

# ============================================================
# Step 4: Summary
# ============================================================
echo ""
echo "===== $(date) Scaffolding Summary ====="

FINAL_SCAFFOLDS="${OUTDIR}/${PREFIX}_scaffolds_final.fa"

if [ -f "${FINAL_SCAFFOLDS}" ]; then
    awk '/^>/{if(n++) print len; len=0; next} {len+=length($0)} END{print len}' "${FINAL_SCAFFOLDS}" | \
        sort -rn | awk '{vals[NR]=$1; sum+=$1}
        END{
            half=sum/2; cum=0;
            for(i=1;i<=NR;i++){
                cum+=vals[i];
                if(cum>=half && !found){
                    printf "N50: %.2f Mb\n", vals[i]/1e6;
                    found=1
                }
            }
            printf "Longest: %.2f Mb\n", vals[1]/1e6;
            printf "Total: %.2f Gb\n", sum/1e9
        }'
    echo "Scaffolds:  $(grep -c '^>' ${FINAL_SCAFFOLDS})"

    N_TOTAL=$(grep -v "^>" "${FINAL_SCAFFOLDS}" | tr -d 'N' | tr -d '\n' | wc -c)
    TOT=$(grep -v "^>" "${FINAL_SCAFFOLDS}" | tr -d '\n' | wc -c)
    N_COUNT=$((TOT - N_TOTAL))
    GAP_BLOCKS=$(grep -v "^>" "${FINAL_SCAFFOLDS}" | tr -d '\n' | grep -oP 'N+' | wc -l)
    echo "N bases:     ${N_COUNT}"
    echo "Gap blocks:  ${GAP_BLOCKS}"

    echo ""
    echo "--- Scaffolds >= 10 Mb ---"
    awk '/^>/{if(n++) print len; len=0; next} {len+=length($0)} END{print len}' "${FINAL_SCAFFOLDS}" | \
        sort -rn | awk '{if($1>=10000000){c++; sum+=$1}} END{printf "%d scaffolds >= 10 Mb = %.2f Gb\n", c, sum/1e9}'

    echo ""
    echo "Output: ${FINAL_SCAFFOLDS}"
else
    echo "[$(date)] ERROR: Scaffold output not found!"
    exit 1
fi

echo ""
echo "===== $(date) YaHS v3 COMPLETE ====="
