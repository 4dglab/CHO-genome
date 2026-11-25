##################################
#!/bin/python3
# 2026-08-21
# CHO_ONT HiC-Pro
# env: hic
# REQUIRES: source ../../config.sh before running snakemake
###################################

import os

SAMPLE = "K1_merged"
TRIM = os.environ["HICPRO_TRIM_MK"]
CONFIG = os.environ["HICPRO_CONFIG"]
MAT_FORMAT = os.environ["HICPRO_MAT_MK"]
SIZE = os.environ["GENOME_SIZES"]
RAW = os.environ["HIC_RAW_DIR"]
CLEAN = os.environ["HIC_CLEAN_DIR"]

rule all:
    input:
        "hic-pro.done",
        expand("mat_format_trans/CHO_ONT/{sample}.mcool", sample=SAMPLE)

rule bl_trimming:
    input:
        f"{RAW}/{{sample}}/{{sample}}_R1.fq.gz",
        f"{RAW}/{{sample}}/{{sample}}_R2.fq.gz"
    output:
        f"{CLEAN}/{{sample}}/{{sample}}_R1.fq.gz",
        f"{CLEAN}/{{sample}}/{{sample}}_R2.fq.gz"
    params:
        sample_dir = f"{RAW}/{{sample}}",
        out_dir = f"{CLEAN}/{{sample}}"
    shell:
        """
        make -f {TRIM} filter_trim sample_dir={params.sample_dir} outdir={params.out_dir}
        rm -f {params.out_dir}/*_clean_R*.fq.gz.*
        """

rule change_name:
    input:
        r1 = f"{CLEAN}/{{sample}}/{{sample}}_R1.fq.gz",
        r2 = f"{CLEAN}/{{sample}}/{{sample}}_R2.fq.gz"
    output:
        r1 = f"{CLEAN}/{{sample}}/{{sample}}_R1.fastq.gz",
        r2 = f"{CLEAN}/{{sample}}/{{sample}}_R2.fastq.gz"
    shell:
        """
        mv {input.r1} {output.r1}
        mv {input.r2} {output.r2}
        """

rule hic_pro:
    input:
        expand(f"{CLEAN}/{{sample}}/{{sample}}_R1.fastq.gz", sample=SAMPLE)
    output:
        touch("hic-pro.done")
    params:
        indir = f"{CLEAN}",
        outdir = "hic-pro_out"
    shell:
        """
        if [ -d "{params.outdir}" ]; then
            mv "{params.outdir}" "{params.outdir}.backup.$(date +%Y%m%d%H%M%S)"
        fi
        HiC-Pro -i "{params.indir}" -o "{params.outdir}" -c "{CONFIG}"
        touch {output}
        """

rule valid2mcool:
    input:
        "hic-pro.done"
    output:
        "mat_format_trans/CHO_ONT/{sample}.mcool"
    params:
        "hic-pro_out/hic_results/data/{sample}/{sample}.allValidPairs",
        "mat_format_trans/CHO_ONT/{sample}"
    shell:
        "make -f {MAT_FORMAT} \
        chromsizes={SIZE} \
        op={params[1]} \
        valid={params[0]} \
        valid2cool"
