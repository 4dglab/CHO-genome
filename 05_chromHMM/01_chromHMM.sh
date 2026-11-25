#!/usr/bin/env bash
#####################################
# 2026-08-31
# chromHMM for chip-CriGri_PUMC
# env: chromhmm-env
# REQUIRES (set in ../config.sh):
#   GENOME_FA             reference genome fasta
#   GENOME_SIZES          chrom sizes file
#   GENOME_GTF            gene annotation (a genePred is generated from it)
#   CHIP_BAM_DIR          sorted ChIP-seq BAMs
#   CHROMHMM_PROJECT_DIR  chromHMM working directory
#   CHROMHMM_HOME         ChromHMM install (ChromHMM.jar, CHROMSIZES/, COORDS/)
#   CHROMHMM_GENOME       genome name used inside ChromHMM
#   GTF_TO_GENEPRED       gtfToGenePred executable
#   CHROMHMM_JAVA_OPTS    Java heap and GC options
#####################################

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"

project_dir="${CHROMHMM_PROJECT_DIR}"
genome_file="${GENOME_FA}"
chrom_size="${GENOME_SIZES}"
sorted_bam_dir="${CHIP_BAM_DIR}"
input_bam_dir="${project_dir}/chromHMM_inputBAM"
chromHMM_file="${CHROMHMM_HOME}/ChromHMM.jar"
gtf_file="${GENOME_GTF}"
genePred_file="${project_dir}/$(basename "${GENOME_GTF%.*}").genePred"

genome_name="${CHROMHMM_GENOME}"

# set -euo pipefail

# # 1. Prepare the input file structure
mkdir -p $input_bam_dir
# mv -v $sorted_bam_dir/*.bam $input_bam_dir
cd $project_dir

ln -s $sorted_bam_dir/*.bam \
  $input_bam_dir


# 2. BinarizeBam input files
# All paths in the configuration TSV must be relative; otherwise Java will
# concatenate them incorrectly.
echo "doing binarizeBam..."
java ${CHROMHMM_JAVA_OPTS} \
  -jar $chromHMM_file BinarizeBam \
  -b 200 -f 2 -t chromHMM_sig_out \
  "${CHROMHMM_HOME}/CHROMSIZES/${genome_name}.txt" \
  chromHMM_inputBAM \
  mark_table.tsv \
  chromHMM_out/count


# # 3. Build the annotation file for a non-model species
"${GTF_TO_GENEPRED}" -genePredExt "${gtf_file}" "${genePred_file}"

coord_root="${CHROMHMM_HOME}/COORDS"
mkdir -p "$coord_root/$genome_name"

java -jar $chromHMM_file ConvertGeneTable -nobin -noheader \
  -u "$coord_root" \
  $genePred_file $genome_name $genome_name


# 4. Build a chromatin state model (10-state model)
echo "doing LearnModel..."
java ${CHROMHMM_JAVA_OPTS} \
  -jar $chromHMM_file LearnModel \
  -noimage \
  -b 200 -l $chrom_size -p 20 \
  chromHMM_out/count chromHMM_out/model 10 $genome_name

# 5. Analyze overlap with known annotations
java ${CHROMHMM_JAVA_OPTS} \
  -jar $chromHMM_file OverlapEnrichment \
  -b 200 chromHMM_out/model/${genome_name}_10_segments.bed \
  "${coord_root}/${genome_name}" \
  chromHMM_out/enrichment/enrich