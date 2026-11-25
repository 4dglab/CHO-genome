#!/usr/bin/env bash
# env: hic
# 2026-08-20
# REQUIRES (set in ../../config.sh):
#   HIC_PREFIX        genome prefix used by the HiC-Pro run
#   HIC_REF_DIR       directory for the clean genome, sizes and restriction BED
#   BT2_INDEX_DIR     bowtie2 index root
#   HICPRO_ANNO_DIR   step1_all_chr_size.pl and pureLength_fasta.py
#   HICPRO_HOME       HiC-Pro install (bin/utils/digest_genome.py)
# ==============!!input contents below!!===================
# input:
# genome name
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

prefix="${HIC_PREFIX}"
prefix_clean="${prefix}_clean"
# dir to save ref files
ref_dir="${HIC_REF_DIR}"
bt2_dir="${BT2_INDEX_DIR}"
# no-chr-segments' prefix
# e.g."ptg" "GL|JH|MU|chrM"
segments_prefix=""
# restriction enzyme
restriction_enzyme="HaeIII"
restriction_sites=GG^CC

max_thread_num=20
# ===============!!do not change contents below!!=======================
# files
raw_genome_file="$ref_dir/${prefix}.fa"
all_chrom_sizes_file="$ref_dir/${prefix}_all_chrom_sizes"
clean_genome_file="$ref_dir/${prefix_clean}.fa"
restriction_out_file="$ref_dir/${prefix_clean}_${restriction_enzyme}_frag.bed"



# tool scripts
sizes_pl="${HICPRO_ANNO_DIR}/step1_all_chr_size.pl"
filter_script="${HICPRO_ANNO_DIR}/pureLength_fasta.py"
restriction_py="${HICPRO_HOME}/bin/utils/digest_genome.py"



# 0.all chr sizes stat
# perl $sizes_pl $raw_genome_file $all_chrom_sizes_file # all_chrom_sizes_file is the output!! difined by youself!!
echo "==============="
echo "raw .fa chrom sizes stat finished."
echo "==============="

# # 1. get clean .fa
# ordinary:
# python $filter_script -g $raw_genome_file \
#   -p $prefix_clean -a $segments_prefix -o $ref_dir

# in cho_ont:
samtools faidx "$raw_genome_file"

for i in {1..11}; do
    samtools faidx "$raw_genome_file" "scaffold_$i"
done > "$clean_genome_file"

echo "clean .fa is ${ref_dir}/${prefix_clean}.fa"
echo "==============="

# 2. stat clean chr sizes
perl $sizes_pl $clean_genome_file "$ref_dir/${prefix_clean}.sizes"
echo "clean .fa chrom sizes stat finished."
echo "==============="

# 3. bowtie2-build
bowtie2-build --threads $max_thread_num $clean_genome_file "${bt2_dir}/${prefix_clean}"
echo "bt2 index finished at ${bt2_dir}/${prefix_clean}"
echo "==============="

# 4. restriction_sites
python3 $restriction_py $clean_genome_file -r $restriction_sites \
  -o $restriction_out_file
echo "Genereated restriction file: ${restriction_out_file}"
echo "==============="