#!/bin/bash
# env: RNA-seq
# 2026-08-28
# REQUIRES (set in ../../../config.sh):
#   ANNOT_REF_DIR      reference root: <genome>/<genome>.fa and <genome>.gtf
#   HISAT2_INDEX_DIR   output directory for the HISAT2 indices
#   GENOME_REF_DIR     CriGri_PUMC.fa and CriGri_PUMC.gtf for the last index

# gffread ${ANNOT_REF_DIR}/cg13b/ncbi_dataset/data/cg13b/cg13b.gff \
# 	-T -E -o ${ANNOT_REF_DIR}/cg13b/ncbi_dataset/data/cg13b/cg13b.gtf

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../../ && pwd)/config.sh"

# known genomes
GENOME=("cg13a" "cg13b" "cg20")
BASE="${ANNOT_REF_DIR}"
HT2_DIR="${HISAT2_INDEX_DIR}"
THREADS=16

for genome in "${GENOME[@]}"
do
    echo "processing ${genome}"

    GTF="${BASE}/${genome}/${genome}.gtf"
    FA="${BASE}/${genome}/${genome}.fa"
    SS="${BASE}/${genome}/${genome}.ss"
    EXON="${BASE}/${genome}/${genome}.exon"
    INDEX="${HT2_DIR}/${genome}"

    hisat2_extract_splice_sites.py "$GTF" > "$SS"
    hisat2_extract_exons.py "$GTF" > "$EXON"

    hisat2-build -p "$THREADS" \
        --ss "$SS" \
        --exon "$EXON" \
        "$FA" \
        "$INDEX"
    echo "index for ${genome} is finished."
done

# CriGri_PUMC
GENOME_BASE="${GENOME_REF_DIR}"
PUMC_FA="${GENOME_BASE}/CriGri_PUMC.fa"
PUMC_INDEX="${HISAT2_INDEX_DIR}/CriGri_PUMC"

hisat2_extract_splice_sites.py \
    "${GENOME_BASE}/CriGri_PUMC.gtf" > "${GENOME_BASE}/CriGri_PUMC.ss"
 hisat2_extract_exons.py \
    "${GENOME_BASE}/CriGri_PUMC.gtf" > "${GENOME_BASE}/CriGri_PUMC.exon"

hisat2-build -p "$THREADS" \
    --ss "${GENOME_BASE}/CriGri_PUMC.ss" \
    --exon "${GENOME_BASE}/CriGri_PUMC.exon" \
    "$PUMC_FA" \
    "$PUMC_INDEX"
echo "index for CriGri_PUMC is finished."


