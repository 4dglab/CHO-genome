#!/bin/bash
# env: RNA-seq
# 2026-08-28
# REQUIRES (set in ../../../config.sh):
#   ANNOT_REF_DIR      reference root: <genome>/<genome>.fa and <genome>.gtf
#   HISAT2_INDEX_DIR   output directory for the HISAT2 indices
#   ONT_REF_DIR        cho_ont.fa and cho_ont.gtf for the last index

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

# cho_ont
ONT_BASE="${ONT_REF_DIR}"
ONT_FA="${ONT_BASE}/cho_ont.fa"
ONT_INDEX="${HISAT2_INDEX_DIR}/cho_ont"

hisat2_extract_splice_sites.py \
    "${ONT_BASE}/cho_ont.gtf" > "${ONT_BASE}/cho_ont.ss"
 hisat2_extract_exons.py \
    "${ONT_BASE}/cho_ont.gtf" > "${ONT_BASE}/cho_ont.exon"

hisat2-build -p "$THREADS" \
    --ss "${ONT_BASE}/cho_ont.ss" \
    --exon "${ONT_BASE}/cho_ont.exon" \
    "$ONT_FA" \
    "$ONT_INDEX"
echo "index for cho_ont is finished."


