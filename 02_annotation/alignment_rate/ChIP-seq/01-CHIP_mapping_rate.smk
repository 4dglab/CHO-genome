#!/usr/bin/env python3
# 2026-08-22
# env: CHIP-seq
# REQUIRES: source ../../../config.sh before running snakemake

import os

# cg20 stands for CriGri-PICRH-1.0 (2020)
GENOME = ["cho_ont", "cg_20"]
CATAGORY = "CHO_2021"
SAMPLE = ["H3K4me3","H3K27ac","H3K9me3","input"]
TIME = ["T3","T7"]
data_dir = os.environ["CHIP_ALIGN_DIR"]
BT2_INDEX = os.environ["BT2_INDEX_DIR"]

wildcard_constraints:
    catagory="CHO_2021",
    sample="(?:H3K4me3|H3K27ac|H3K9me3|input)",
    time="(?:T3|T7)",
    genome="(?:cho_ont|cg_20)"

rule all:
    input:
        # expand("bam/{catagory}/{sample}_{time}_{genome}_sorted.bam.csi", sample=SAMPLE, time=TIME, catagory=CATAGORY),
        # expand("macs3_result/{catagory}/{sample}/{sample}_{time}_macs3_peaks.narrowPeak", sample=SAMPLE, time=TIME, catagory=CATAGORY),
        expand("bam/{catagory}/{sample}_{time}_{genome}.sam", sample=SAMPLE, time=TIME, catagory=CATAGORY,genome=GENOME)

rule fastp:
    input:
        f"{data_dir}/{{catagory}}/{{sample}}_{{time}}.fastq.gz"
    output:
        "clean_data/{catagory}/{sample}_{time}_clean.fastq.gz"
    log:
        "clean_data/{catagory}/{sample}_{time}_fastp.log"
    shell:
        "fastp -i {input} -o {output} > {log} 2>&1"

rule bowtie2_mapping:
    input:
        "clean_data/{catagory}/{sample}_{time}_clean.fastq.gz"
    output:
        "bam/{catagory}/{sample}_{time}_{genome}.sam"
    log:
        "bam/{catagory}/{sample}_{time}_{genome}_bowtie2.log"
    params:
        f"{BT2_INDEX}/{{genome}}/{{genome}}"
    shell:
        "bowtie2 -p 10 -x {params} -U {input} -S {output} > {log} 2>&1"
