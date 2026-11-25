#!/usr/bin/env bash
# PURPOSE: NextPolish2 first-round polishing
# INPUT:   hifiasm p_ctg + HiFi + NGS clean + yak k21/k31
# OUTPUT:  cho_hifiasm_v3_polished1.fa
# TOOL:    NextPolish2 0.2.2 / meryl / winnowmap
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -eo pipefail


# ============================================================
# NextPolish2 v3: First polish on raw hifiasm v3 contigs
# ============================================================
# Plan step: hifiasm → Polish_1 → purge_dups → YaHS → ...
#
# Reuses:
#   - NGS clean reads (from nextpolish2/ngs_clean/)
#   - yak k-mer DBs (nextpolish2/kmer/k21.yak, k31.yak) — same NGS reads
#   - HiFi reads
# Builds fresh:
#   - meryl count + print on raw contigs
#   - winnowmap HiFi → raw contigs
#   - nextPolish2 polish
# ============================================================

WORKDIR="${WORK_DIR}/nextpolish2_v3"
GENOME="${WORK_DIR}/hifiasm/v3_diploid/CH6F_G1_hifiasm.bp.p_ctg_wrapped.fasta"
HIFI1="${WORK_DIR}/hifi_fastq/CH6F_G1-1_hifi.fq.gz"
HIFI2="${WORK_DIR}/hifi_fastq/CH6F_G1-2_hifi.fq.gz"
NGS_R1="${WORK_DIR}/nextpolish2/ngs_clean/R1.clean.fq.gz"
NGS_R2="${WORK_DIR}/nextpolish2/ngs_clean/R2.clean.fq.gz"
YAK21="${WORK_DIR}/nextpolish2/kmer/k21.yak"
YAK31="${WORK_DIR}/nextpolish2/kmer/k31.yak"

THREADS=48

cd "${WORKDIR}"

LOG="${WORKDIR}/logs/run.log"
exec > >(tee -a "${LOG}") 2>&1

echo "========== $(date) NextPolish2 v3 (Polish_1) Start =========="
echo "Input genome: ${GENOME}"
echo "HiFi reads:   ${HIFI1}"
echo "              ${HIFI2}"
echo "Threads:      ${THREADS}"

# Verify inputs
for f in "${GENOME}" "${HIFI1}" "${HIFI2}" "${NGS_R1}" "${NGS_R2}" "${YAK21}" "${YAK31}"; do
    if [ ! -f "$f" ]; then
        echo "[$(date)] ERROR: Missing input: $f"
        exit 1
    fi
done
echo "[$(date)] All input files verified"

# Quick stats on input contigs
N_CTG=$(grep -c "^>" "${GENOME}")
TOTAL_BP=$(grep -v "^>" "${GENOME}" | tr -d '\n' | wc -c)
echo "[$(date)] Input: ${N_CTG} contigs, ${TOTAL_BP} bp"

# ============================================================
# Step 1: meryl count k=15 on raw contigs
# ============================================================
echo "[$(date)] Step 1: meryl count k=15 on raw contigs..."

conda activate nextpolish2-env

rm -rf "${WORKDIR}/hifi_map/merylDB"
meryl count k=15 \
    output "${WORKDIR}/hifi_map/merylDB" \
    "${GENOME}" \
    threads=${THREADS} \
    2>&1 | tee "${WORKDIR}/logs/meryl_count.log"

echo "[$(date)] meryl count DONE"

# ============================================================
# Step 2: meryl print repetitive k-mers
# ============================================================
echo "[$(date)] Step 2: meryl print repetitive k-mers..."

meryl print greater-than distinct=0.9998 \
    "${WORKDIR}/hifi_map/merylDB" \
    > "${WORKDIR}/hifi_map/repetitive_k15.txt" \
    2> "${WORKDIR}/logs/meryl_print.log"

N_REP=$(wc -l < "${WORKDIR}/hifi_map/repetitive_k15.txt")
echo "[$(date)] meryl print DONE: ${N_REP} repetitive k-mers"

# ============================================================
# Step 3: winnowmap — HiFi reads → raw contigs
# ============================================================
echo "[$(date)] Step 3: winnowmap — mapping HiFi reads to raw contigs..."

winnowmap \
    -t ${THREADS} \
    -W "${WORKDIR}/hifi_map/repetitive_k15.txt" \
    -ax map-pb \
    "${GENOME}" \
    "${HIFI1}" "${HIFI2}" \
    2> "${WORKDIR}/logs/winnowmap.log" \
    | samtools sort -@ 16 -o "${WORKDIR}/hifi_map/hifi.map.sort.bam" -

echo "[$(date)] winnowmap DONE: $(du -h ${WORKDIR}/hifi_map/hifi.map.sort.bam | cut -f1)"

# ============================================================
# Step 4: samtools index (CSI)
# ============================================================
echo "[$(date)] Step 4: samtools index..."

rm -rf "${WORKDIR}/hifi_map/merylDB"

samtools index -c "${WORKDIR}/hifi_map/hifi.map.sort.bam"

echo "[$(date)] samtools index DONE"

# ============================================================
# Step 5: nextPolish2
# ============================================================
echo "[$(date)] Step 5: nextPolish2 polishing..."

nextPolish2 \
    -t ${THREADS} \
    "${WORKDIR}/hifi_map/hifi.map.sort.bam" \
    "${GENOME}" \
    "${YAK21}" \
    "${YAK31}" \
    > "${WORKDIR}/output/cho_hifiasm_v3_polished1.fa" \
    2> "${WORKDIR}/logs/nextpolish2.log"

echo "[$(date)] nextPolish2 DONE"

# ============================================================
# Step 6: QC
# ============================================================
echo "[$(date)] Step 6: QC..."

FA="${WORKDIR}/output/cho_hifiasm_v3_polished1.fa"

N_CTG_OUT=$(grep -c "^>" "${FA}")
TOTAL_OUT=$(grep -v "^>" "${FA}" | tr -d '\n' | wc -c)
N_BASES=$(grep -v "^>" "${FA}" | tr -d 'N' | tr -d '\n' | wc -c)
N_COUNT=$((TOTAL_OUT - N_BASES))

echo ""
echo "========== $(date) NextPolish2 v3 COMPLETE =========="
echo "Input:          ${N_CTG} contigs, ${TOTAL_BP} bp"
echo "Polished:       ${N_CTG_OUT} contigs, ${TOTAL_OUT} bp"
echo "N bases:        ${N_COUNT}"
echo "Output:         ${FA}"
