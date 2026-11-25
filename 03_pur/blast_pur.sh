#!/usr/bin/env bash
# PURPOSE: PUR cross-species BLAST (mouse/rat)
# INPUT:   PUR sequences
# OUTPUT:  BLAST results
# TOOL:    BLAST+ 2.17.0
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"


# ============================================================
# ============================================================
set -euo pipefail


WORK="${WORK_DIR}/reviewer_response/pur_analysis"
PDIR="${WORK}/pur_full_genome"
BDIR="${WORK}/blast_pur"

QRY="${WORK_DIR}/nextpolish2_review4_gapfilled/block46_translocation/output/cho_review5_gapfilled_polished.fa"
QRY_FAI="${QRY}.fai"
MOUSE_RAW="${RAW_DIR}/genome_assembly/assemblyref/GCF_000001635.27_GRCm39_genomic.fna"
RAT_RAW="${RAW_DIR}/genome_assembly/assemblyref/GRCr8/GCF_036323735.1_GRCr8_genomic.fna"
THREADS=16
EVALUE=1e-5
OUTFMT="6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send"

mkdir -p "${BDIR}"
cd "${BDIR}"

echo "===== $(date) Phase C: cross-species homology ====="

PUR_FA="${PDIR}/pur_full.ge1kb.fa"
if [ ! -s "${PUR_FA}" ]; then
    echo "[$(date)] extracting PUR >=1kb..."
    python3 "${WORK}/scripts/04_extract_pur_ge1kb.py" \
        --region "${PDIR}/pur_full.region" \
        --fasta "${QRY}" --fai "${QRY_FAI}" \
        --outfa "${PUR_FA}" \
        --outbed "${PDIR}/pur_full.ge1kb.bed"
fi
echo "[$(date)] PUR >=1kb: $(grep -c '^>' ${PUR_FA}) seqs"

conda activate busco-env

MOUSE_FA="${BDIR}/GRCm39_genomic.fna"
if [ -f "${MOUSE_RAW}.gz" ]; then
    gunzip -c "${MOUSE_RAW}.gz" > "${MOUSE_FA}"
elif [ -f "${MOUSE_RAW}" ]; then
    cp "${MOUSE_RAW}" "${MOUSE_FA}"
fi

if [ -s "${MOUSE_FA}" ]; then
    SIZE=$(stat -c %s "${MOUSE_FA}")
    if [ "${SIZE}" -lt 500000000 ]; then
        echo "[WARN] GRCm39 file suspiciously small (${SIZE} bytes) - mouse blastn SKIPPED. Please re-download."
    else
        echo "[$(date)] makeblastdb mouse..."
        makeblastdb -in "${MOUSE_FA}" -dbtype nucl -out db_mouse_GRCm39
        echo "[$(date)] blastn mouse..."
        blastn -task megablast -query "${PUR_FA}" -db db_mouse_GRCm39 \
            -out pur.vs.mouse.tsv -outfmt "${OUTFMT}" \
            -evalue ${EVALUE} -max_target_seqs 1 -num_threads ${THREADS}
        echo "[$(date)] mouse hits: $(wc -l < pur.vs.mouse.tsv)"
    fi
else
    echo "[WARN] GRCm39 missing - mouse blastn SKIPPED."
fi

RAT_FA="${BDIR}/GRCr8_genomic.fna"
if [ -f "${RAT_RAW}.gz" ]; then
    gunzip -c "${RAT_RAW}.gz" > "${RAT_FA}"
elif [ -f "${RAT_RAW}" ]; then
    cp "${RAT_RAW}" "${RAT_FA}"
fi

if [ -s "${RAT_FA}" ]; then
    echo "[$(date)] makeblastdb rat..."
    makeblastdb -in "${RAT_FA}" -dbtype nucl -out db_rat_GRCr8
    echo "[$(date)] blastn rat..."
    blastn -task megablast -query "${PUR_FA}" -db db_rat_GRCr8 \
        -out pur.vs.rat.tsv -outfmt "${OUTFMT}" \
        -evalue ${EVALUE} -max_target_seqs 1 -num_threads ${THREADS}
    echo "[$(date)] rat hits: $(wc -l < pur.vs.rat.tsv)"
else
    echo "[WARN] GRCr8 missing - rat blastn SKIPPED."
fi

echo "===== $(date) Phase C blast COMPLETE ====="
