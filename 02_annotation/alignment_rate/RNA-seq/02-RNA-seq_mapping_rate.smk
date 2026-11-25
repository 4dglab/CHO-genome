#!/usr/bin/env python3
# 2026-08-22
# env: RNA-seq
# REQUIRES: source ../../../config.sh before running snakemake

import os

# fastp quality control and adapter trimming
# cg13a stands for Cgr1.0 (2013), cg13b stands for C_griseus_v1.0 (2013), cg20 stands for CriGri-PICRH-1.0 (2020), cho_ont stands for our genome assembly (2026)
SAMPLE = ["sdu_1","sdu_2","sdu_3"]
GENOME = ["cg13a","cg13b","cg20","cho_ont"]
data_dir = os.environ["RNA_ALIGN_DIR"]
HISAT2_INDEX = os.environ["HISAT2_INDEX_DIR"]

rule all:
    input:
        # expand("fastp/{sample}_1.fastq.gz", sample=config["samples"]),
        # expand("fastp/{sample}_2.fastq.gz", sample=config["samples"])
        expand("mapping/{sample}/{sample}_{genome}.sam", sample=SAMPLE, genome=GENOME)
        # expand("Mapping/{sample}/{sample}_{genome}_sorted.bam.bai", sample=SAMPLE, genome=GENOME),
        # expand("cufflinks_res/{sample}/genes.fpkm_tracking", sample=SAMPLE),
        # expand("Mapping/{sample}/{sample}_{genome}.mapping_chr_stats", sample=SAMPLE, genome=GENOME)


rule fastp:
    input:
        f"{data_dir}/{{sample}}/{{sample}}_R1.fastq.gz",
        f"{data_dir}/{{sample}}/{{sample}}_R2.fastq.gz"
    output:
        "clean_data/{sample}/{sample}_R1_clean.fastq.gz",
        "clean_data/{sample}/{sample}_R2_clean.fastq.gz"
    log:
        "clean_data/{sample}/{sample}.log"
    shell:
        "fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w 10 > {log} 2>&1"

rule ht2_mapping:
    input:
        "clean_data/{sample}/{sample}_R1_clean.fastq.gz",
        "clean_data/{sample}/{sample}_R2_clean.fastq.gz"
    output:
        "mapping/{sample}/{sample}_{genome}.sam"
    params:
        f"{HISAT2_INDEX}/{{genome}}/{{genome}}"
    log:
        "mapping/{sample}/{sample}_{genome}.log"
    shell:
        "hisat2 -1 {input[0]} -2 {input[1]} -p 10 -x {params} -S {output} > {log} 2>&1"

# rule samtools_sort:
#     input:
#         "Mapping/{sample}/{sample}_{genome}.sam"
#     output:
#         "Mapping/{sample}/{sample}_{genome}_sorted.bam"
#     log:
#         "Mapping/{sample}/{sample}_{genome}_sort.log"
#     params:
#         "Mapping/{sample}/{sample}_{genome}_sorted.bam.temp"
#     shell:
#         "samtools sort -O BAM -@ 10 -m 2G -o {output} -T {params} {input} > {log} 2>&1"

# View mapping statistics for each chromosome
# rule mapping_stats:
#     input:
#         "Mapping/{sample}/{sample}_{genome}_sorted.bam"
#     output:
#         "Mapping/{sample}/{sample}_{genome}.mapping_chr_stats"
#     shell:
#         "samtools idxstats -@ 10 {input} > {output}"

#  rule bam_index:
# 	input:
# 		"Mapping/{sample}/{sample}_{genome}_sorted.bam"
# 	output:
# 		"Mapping/{sample}/{sample}_{genome}_sorted.bam.bai"
# 	shell:
# 		"samtools index {input}"

# rule cufflinks:
# 	input:
# 		"Mapping/{sample}/{sample}_{genome}_sorted.bam"
# 	output:
# 		"cufflinks_res/{sample}/genes.fpkm_tracking"
# 	params:
# 		"cufflinks_res/{sample}"
# 	log:
# 		"cufflinks_res/{sample}/cufflinks.log"
# 	shell:
# 		"cufflinks -o {params} -p 10 -G {gtf} {input} > {log} 2>&1"
