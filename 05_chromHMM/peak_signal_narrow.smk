# Call narrow peaks.
# env: CHIP-seq
# REQUIRES: source ../config.sh before running snakemake

import os

SAMPLE = ["H3K4me3","H3K27ac"]
BAM_DIR = os.environ["CHIP_BAM_DIR"]
BW_DIR = os.environ["HISTONE_BW_DIR"]
LOG_DIR = os.environ["MACS3_OUT_DIR"]
chrom_sizes = os.environ["GENOME_SIZES"]

rule all:
    input:
        # expand("macs3_out/{sample}/{sample}_peaks.narrowPeak",sample=SAMPLE)
        expand(f"{BW_DIR}/{{sample}}/{{sample}}_FE.bw",sample=SAMPLE)

rule macs3_callpeak:
    input:
        f"{BAM_DIR}/{{sample}}_sorted.bam"
    output:
        f"{BW_DIR}/{{sample}}/{{sample}}_treat_pileup.bdg",
        f"{BW_DIR}/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        f"{BW_DIR}/{{sample}}/",
        f"{BAM_DIR}/input_sorted.bam"
    log:
        f"{LOG_DIR}/{{sample}}.log"
    shell:
        "macs3 callpeak \
        -c {params[1]} -t {input} \
        -f BAM -g 2437314528 \
        -n {wildcards.sample} --outdir {params[0]} \
        -q 0.05 -m 3 50 -B --SPMR > {log} 2>&1"

rule foldenrichment:
    input:
        f"{BW_DIR}/{{sample}}/{{sample}}_treat_pileup.bdg",
        f"{BW_DIR}/{{sample}}/{{sample}}_control_lambda.bdg"
    output:
        f"{BW_DIR}/{{sample}}/{{sample}}_FE.bdg"
    shell:
        "macs3 bdgcmp \
        -t {input[0]} -c {input[1]} -m FE -o {output}"

# rule bedtools_slop:
#     input:
#         fe_bdg = f"{BW_DIR}/{{sample}}/{{sample}}_FE.bdg",
#         chrsz  = chrom_sizes
#     output:
#         cleaned = "{sample}.fc.clean.bedgraph"
#     shell:
#         r"""
#         bedtools slop -i {input.fe_bdg} -g {input.chrsz} -b 0 \
#         | awk '{{ if ($3 != -1) print $0 }}' \
#         | bedClip stdin {input.chrsz} {output.cleaned}
#         """
rule sortbdg:
    input:
         f"{BW_DIR}/{{sample}}/{{sample}}_FE.bdg"
    output:
         f"{BW_DIR}/{{sample}}/{{sample}}_FE_sorted.bdg"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule bdg2bw:
    input:
        f"{BW_DIR}/{{sample}}/{{sample}}_FE_sorted.bdg"
    output:
        f"{BW_DIR}/{{sample}}/{{sample}}_FE.bw"
    shell:
        "bedGraphToBigWig \
        {input} {chrom_sizes} {output}"
