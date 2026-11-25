#!/usr/bin/env bash
# PURPOSE: purge_dups haplotypic redundancy removal (439 -> 119 contigs)
# INPUT:   cho_hifiasm_v3_polished1.fa + HiFi
# OUTPUT:  purged_clean.fa
# TOOL:    purge_dups / minimap2
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


set -eo pipefail

conda activate purge_dups-env

# ============================================================
# purge_dups v3 — on polish_1 contig (439 contigs, 0 gaps)
# ============================================================
# Plan step: hifiasm → Polish_1 → purge_dups → YaHS → ...
# ref: plan.md, row 4

WORKDIR="${WORK_DIR}/purge_dups_v3"
GENOME="${WORK_DIR}/nextpolish2_v3/output/cho_hifiasm_v3_polished1.fa"
HIFI_1="${WORK_DIR}/hifi_fastq/CH6F_G1-1_hifi.fq.gz"
HIFI_2="${WORK_DIR}/hifi_fastq/CH6F_G1-2_hifi.fq.gz"
THREADS=48
PREFIX="cho_polished1"

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

LOG="${WORKDIR}/logs/run.log"
exec > >(tee -a "${LOG}") 2>&1

echo ">>>> $(date) purge_dups v3 workflow start <<<<"
echo "Target:   ${GENOME}"
echo "Reads:    ${HIFI_1}"
echo "          ${HIFI_2}"
echo "Threads:  ${THREADS}"

# Verify inputs
for f in "${GENOME}" "${HIFI_1}" "${HIFI_2}"; do
    if [ ! -f "$f" ]; then
        echo "[$(date)] ERROR: Missing input: $f"
        exit 1
    fi
done
echo "[$(date)] All inputs verified"

N_CTG=$(grep -c "^>" "${GENOME}")
TOTAL_BP=$(grep -v "^>" "${GENOME}" | tr -d '\n' | wc -c)
echo "[$(date)] Input: ${N_CTG} contigs, ${TOTAL_BP} bp"

# ============================================================
# Step 1: minimap2 — HiFi reads → contigs → PAF
# ============================================================
echo "[$(date)] Step 1: minimap2 mapping HiFi reads..."

minimap2 \
    -x map-hifi \
    -t ${THREADS} \
    "${GENOME}" \
    "${HIFI_1}" "${HIFI_2}" \
    > "${PREFIX}.paf" \
    2> logs/minimap2_reads.log

echo "[$(date)] minimap2 reads PAF done: $(du -h ${PREFIX}.paf | cut -f1)"

# ============================================================
# Step 2: pbcstat
# ============================================================
echo "[$(date)] Step 2: pbcstat..."

grep -v "^#" "${PREFIX}.paf" > "${PREFIX}.clean.paf"
echo "  Cleaned PAF: $(du -h ${PREFIX}.clean.paf | cut -f1)"

pbcstat "${PREFIX}.clean.paf" \
    2> logs/pbcstat.log

echo "[$(date)] pbcstat done"
ls -lh PB.stat PB.base.cov

# ============================================================
# Step 3: calcuts -d 1 (haploid mode)
# ============================================================
echo "[$(date)] Step 3: calcuts..."

calcuts -d 1 PB.stat \
    > cutoffs \
    2> logs/calcuts.log

echo "[$(date)] cutoffs: $(cat cutoffs)"

# ============================================================
# Step 4: split_fa + self-alignment
# ============================================================
echo "[$(date)] Step 4: split_fa + self-alignment..."

split_fa "${GENOME}" \
    > "${PREFIX}.split.fa" \
    2> logs/split_fa.log

N_CHUNKS=$(grep -c '^>' "${PREFIX}.split.fa")
echo "[$(date)] split_fa done: ${N_CHUNKS} chunks"

minimap2 \
    -x asm5 \
    -t ${THREADS} \
    -DP \
    "${PREFIX}.split.fa" \
    "${PREFIX}.split.fa" \
    > "${PREFIX}.split.paf" \
    2> logs/minimap2_self.log

echo "[$(date)] self-alignment done: $(du -h ${PREFIX}.split.paf | cut -f1)"

# ============================================================
# Step 5: purge_dups -2
# ============================================================
echo "[$(date)] Step 5: purge_dups..."

purge_dups \
    -2 \
    -c PB.base.cov \
    -T cutoffs \
    "${PREFIX}.split.paf" \
    > dups.bed \
    2> logs/purge_dups.log

echo "[$(date)] purge_dups done: $(wc -l < dups.bed) duplication entries"

# ============================================================
# Step 6: get_seqs
# ============================================================
echo "[$(date)] Step 6: get_seqs..."

get_seqs \
    dups.bed \
    "${GENOME}" \
    2> logs/get_seqs.log

echo "[$(date)] get_seqs done"

# ============================================================
# Step 7: Summary
# ============================================================
echo ""
echo ">>>> $(date) purge_dups v3 COMPLETE <<<<"
echo ""

ORIG_BASES=$(grep -v "^>" "${GENOME}" | tr -d '\n' | wc -c)
echo "Original:  ${ORIG_BASES} bp"

P_BASES=0
if [ -f purged.fa ]; then
    P_BASES=$(grep -v "^>" purged.fa | tr -d '\n' | wc -c)
    P_CTGS=$(grep -c "^>" purged.fa)
    echo "purged.fa: ${P_CTGS} contigs, ${P_BASES} bp"
else
    echo "purged.fa: NOT GENERATED"
fi

H_BASES=0
if [ -f hap.fa ]; then
    H_BASES=$(grep -v "^>" hap.fa | tr -d '\n' | wc -c)
    H_CTGS=$(grep -c "^>" hap.fa)
    echo "hap.fa:    ${H_CTGS} entries, ${H_BASES} bp  (redundant)"
else
    echo "hap.fa:    NOT GENERATED (no redundancy)"
fi

echo ""
echo "Redundant:  $(awk "BEGIN {printf \"%.2f\", 100*${H_BASES}/${ORIG_BASES}}")%"
