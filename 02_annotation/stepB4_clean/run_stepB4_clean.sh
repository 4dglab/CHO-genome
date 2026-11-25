#!/usr/bin/env bash
# PURPOSE: B4 redundancy removal + pseudogene filtering + evidence tiers L1-L4
# INPUT:   03_agat/genomic.gff
# OUTPUT:  final/genomic.gff/.gtf
# TOOL:    custom
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# Step B4: redundancy removal + pseudogene filtering + evidence classification (L1-L4)
# for the FULL genome (annotation_final). Same logic as the 50Mb test.
#
# Inputs:
#   GENOMIC_GFF  B3 AGAT genomic.gff
#   GENOME_FA    softmasked full genome (for gffread -y)
#   ISOSEQ_BAM   Iso-seq BAM aligned to full genome (from B0c)
#   PSEUDO_LABEL EviAnn pseudo_label GFF
#
# Envs: annotation-env (AGAT/mmseqs), eviann (gffread), HiTE (bedtools),
#       base (matplotlib)
set -euo pipefail

BASE=${WORK_DIR}/annotation_final
STEP=$BASE/stepB4_clean
SCRIPTS=$BASE/scripts
B3=$BASE/stepB_gene/03_agat
EV=$BASE/stepB_gene/02_eviann
PREP=$BASE/stepB_gene/00_prep

GENOMIC_GFF=$B3/genomic.gff
GENOME_FA=$BASE/genome.softmasked.fa
ISOSEQ_BAM=$PREP/isoseq/isoseq_aligned.bam
PSEUDO_LABEL=$EV/genome.fa.pseudo_label.gff
BEDTOOLS=${CONDA_BASE}/envs/HiTE/bin/bedtools

mkdir -p "$STEP"/{logs,primary,cluster,redundancy,pseudogene,evidence,final}
LOGS=$STEP/logs

# Clean previous outputs for idempotent re-run
rm -rf "$STEP"/primary/* "$STEP"/cluster/* "$STEP"/redundancy/* \
       "$STEP"/pseudogene/* "$STEP"/evidence/* "$STEP"/final/*


echo "=== B4-0: remove out-of-bounds models (chunk-boundary artifacts) ==="
python "$SCRIPTS/remove_out_of_bounds.py" \
    --gff "$GENOMIC_GFF" --genome "$GENOME_FA" \
    --out "$STEP/genomic.oob_filtered.gff" \
    --report "$LOGS/out_of_bounds_report.tsv"
GENOMIC_GFF="$STEP/genomic.oob_filtered.gff"
echo "filtered GFF: $(grep -cP '\tgene\t' $GENOMIC_GFF) genes"

echo "=== B4-1: primary representative isoform set (agat keep_longest) ==="
conda activate annotation-env
agat_sp_keep_longest_isoform.pl --cpu 16 -f "$GENOMIC_GFF" \
    -o "$STEP/primary/genomic.primary.gff" >"$LOGS/primary.log" 2>&1
echo "primary: $(grep -cP '\tgene\t' $STEP/primary/genomic.primary.gff) genes"

echo "=== B4-2: extract representative proteins (gffread -y) ==="
conda activate eviann
gffread -y "$STEP/primary/genomic.primary.proteins.fa" \
    -g "$GENOME_FA" "$STEP/primary/genomic.primary.gff" >"$LOGS/gffread.log" 2>&1
echo "proteins: $(grep -c '^>' $STEP/primary/genomic.primary.proteins.fa)"

echo "=== B4-3: mmseqs clustering (>=95% identity, >=80% coverage) ==="
conda activate annotation-env
# mmseqs 16.747c6 segfaults in align step; use newer cached build
MMSEQ=$(find ${CONDA_BASE}/pkgs -maxdepth 3 -type f -name mmseqs -path '*18*' 2>/dev/null | head -1)
[ -z "$MMSEQ" ] && MMSEQ=$(which mmseqs)
echo "mmseqs: $MMSEQ"
"$MMSEQ" easy-cluster "$STEP/primary/genomic.primary.proteins.fa" \
    "$STEP/cluster/clusters" "$STEP/cluster/tmp" \
    --min-seq-id 0.95 -c 0.8 --cov-mode 1 --threads 48 >"$LOGS/mmseqs.log" 2>&1
CLUSTER_TSV="$STEP/cluster/clusters_cluster.tsv"
echo "clusters: $(wc -l < $CLUSTER_TSV) lines"

echo "=== B4-4: redundancy / overlap analysis ==="
python "$SCRIPTS/overlap_analysis.py" \
    --gff "$GENOMIC_GFF" \
    --clusters "$STEP/cluster/clusters_cluster.tsv" \
    --outdir "$STEP/redundancy" 2>&1 | tee "$LOGS/overlap.log"

echo "=== B4-5: pseudogene scan ==="
python "$SCRIPTS/pseudogene_scan.py" \
    --proteins "$STEP/primary/genomic.primary.proteins.fa" \
    --gff "$GENOMIC_GFF" \
    --genome "$GENOME_FA" \
    --pseudo-label "$PSEUDO_LABEL" \
    --outdir "$STEP/pseudogene" 2>&1 | tee "$LOGS/pseudogene.log"

echo "=== B4-6: evidence classification L1-L4 (incl. Iso-seq BAM) ==="
python "$SCRIPTS/evidence_classify.py" \
    --gff "$GENOMIC_GFF" \
    --isoseq-bam "$ISOSEQ_BAM" \
    --bedtools "$BEDTOOLS" \
    --outdir "$STEP/evidence" 2>&1 | tee "$LOGS/evidence.log"

echo "=== B4-7: final GFF + filtering pipeline report ==="
conda activate base
python "$SCRIPTS/filtering_report.py" \
    --gff "$STEP/evidence/genomic.evidence.gff" \
    --pseudo-drop "$STEP/pseudogene/pseudogene_dropped.tsv" \
    --redundant-drop "$STEP/redundancy/dropped_redundant.tsv" \
    --outdir "$STEP/final" 2>&1 | tee "$LOGS/filtering_report.log"

echo "=== B4-8: convert final GFF to GTF ==="
conda activate annotation-env
agat_convert_sp_gff2gtf.pl --cpu 16 --gff "$STEP/final/genomic.gff" \
    -o "$STEP/final/genomic.gtf" >"$LOGS/gff2gtf.log" 2>&1
echo "gtf: $(wc -l < $STEP/final/genomic.gtf) lines"

echo "STEPB4_CLEAN_DONE"
echo "final: $STEP/final/genomic.gff + genomic.gtf"
echo "pipeline table: $STEP/final/filtering_pipeline.tsv"
echo "figure: $STEP/final/filtering_pipeline.png"
