#!/bin/bash
# ==============================================================================
# CHO genome release - global configuration
# This is the only file you need to edit before running the pipeline.
# Set the paths to your downloaded raw data, working directory, and tools below.
# ==============================================================================

# ---- Raw sequencing data (download from GSA: PRJCA028069) ----
# Point each variable to the corresponding downloaded file.
export RAW_DIR="/path/to/raw_data"

# PacBio HiFi reads (BAM)
export HIFI_BAM_1="${RAW_DIR}/CH6F_G1-1.hifi_reads.bam"
export HIFI_BAM_2="${RAW_DIR}/CH6F_G1-2.hifi_reads.bam"

# Oxford Nanopore ultra-long reads (two sequencing runs)
export ONT_UL_1="${RAW_DIR}/PUMC-CHO1.pass.ul.1.fq.gz"
export ONT_UL_2="${RAW_DIR}/PUMC-CHO1.pass.ul.2.fq.gz"

# Illumina Hi-C libraries (two libraries, paired-end)
export HIC_LIB1_R1="${RAW_DIR}/CH6F_G3-1-1_R1.fq.gz"
export HIC_LIB1_R2="${RAW_DIR}/CH6F_G3-1-1_R2.fq.gz"
export HIC_LIB2_R1="${RAW_DIR}/CH6F_G3-1-2_R1.fq.gz"
export HIC_LIB2_R2="${RAW_DIR}/CH6F_G3-1-2_R2.fq.gz"

# Illumina NGS WGS reads
export NGS_R1="${RAW_DIR}/CH6F_G1_R1.fq.gz"
export NGS_R2="${RAW_DIR}/CH6F_G1_R2.fq.gz"

# ---- Working directory (all intermediate and final outputs) ----
export WORK_DIR="/path/to/result"

# ---- Threads ----
export THREADS=48

# ---- Third-party tool paths ----
export JUICER_TOOLS="/path/to/juicer_tools_1.22.01.jar"
export BUSCO_DOWNLOAD="/path/to/busco_downloads"   # BUSCO --offline lineage cache

# ---- conda (used by the Python scripts via os.environ) ----
export CONDA_BASE="/path/to/anaconda3"
export CONDA_PREFIX_GAPFILL="${CONDA_BASE}/envs/gapfill-env"
export CONDA_PREFIX_TGS="${CONDA_BASE}/envs/tgsgapcloser-env"
export CONDA_PREFIX_NEXTPOLISH="${CONDA_BASE}/envs/nextpolish2-env"

# ---- Visualization tools ----
export GENOMESYN="/path/to/GenomeSyn-1.2.7"          # GenomeSyn install directory
export GNUPLOT_BIN="/path/to/gnuplot/bin"            # gnuplot executable dir (mummerplot)
export GNUPLOT_PS_DIR="/path/to/gnuplot/PostScript"  # gnuplot PostScript resources

# ==============================================================================
# Reference genome, annotation and read indices (shared by several modules)
# ==============================================================================
# *_clean files keep scaffold_1 to scaffold_11 only (scaffold_12-74 removed).
# How to build them: 06_Hi-C/README.md, section "Reference preparation".
# If the released assembly and its indices already fit your analysis, point
# these variables at the released files and skip the clean step.
export GENOME_FA="/path/to/ref_genome/CriGri_PUMC_clean.fa"          # requires a .fai index
export GENOME_SIZES="/path/to/ref_genome/CriGri_PUMC_clean.sizes"    # <scaffold> <length>
export GENOME_GTF="/path/to/ref_genome/CriGri_PUMC_clean.gtf"        # same scaffolds as GENOME_FA
export GENOME_GFF="/path/to/genomic.highconfidence.filtered.gff"
export GENOME_BED12="/path/to/Chinese_genome_annotation.bed12"
export GENOME_REF_DIR="/path/to/data/CriGri_PUMC/ref_genome"         # full assembly: CriGri_PUMC.fa + .gtf (flat)
export BT2_INDEX_DIR="/path/to/bowtie2_index"                    # <genome>/<genome>.*
export HISAT2_INDEX_DIR="/path/to/hisat2_index"                  # <genome>/<genome>.*

# ==============================================================================
# 02_annotation — alignment-rate comparison (ChIP-seq / RNA-seq)
# ==============================================================================
export ANNOT_REF_DIR="/path/to/reference"                        # <genome>/<genome>.fa + .gtf
export CHIP_ALIGN_DIR="/path/to/CHIP-seq_alignRate"              # <category>/<sample>_<time>.fastq.gz
export RNA_ALIGN_DIR="/path/to/RNA-seq_alignRate/CHO_sdu"        # <sample>/<sample>_R{1,2}.fastq.gz

# ==============================================================================
# 05_chromHMM — chromatin states
# ==============================================================================
export CHIP_BAM_DIR="/path/to/chrom_status/sorted_bam"           # sorted ChIP-seq BAMs
export CHIP_ORI_BAM_DIR="/path/to/chrom_status/ori_bam"          # original single-end BAMs
export HISTONE_BW_DIR="/path/to/chrom_status/bw"                 # <mark>/<mark>_FE.bw
export MACS3_OUT_DIR="/path/to/chrom_status/macs3_out"
export CHROMHMM_PROJECT_DIR="/path/to/chrom_status/chromHMM_project"
export CHROMHMM_HOME="/path/to/ChromHMM"                         # ChromHMM.jar + CHROMSIZES/ + COORDS/
export CHROMHMM_GENOME="CriGri_PUMC_clean"                           # genome name used inside ChromHMM
export GTF_TO_GENEPRED="/path/to/gtfToGenePred"
export CHROMHMM_JAVA_OPTS="-Xms80g -Xmx80g -XX:ParallelGCThreads=20"

# ==============================================================================
# 06_Hi-C — Hi-C / 3D genome analysis
# ==============================================================================
# Output root for the Hi-C module (loop/, TAD/, AB_compartments/, figures/).
export HI_C_DIR="${WORK_DIR}/hi_c"

# ---- Hi-C input matrices (produced by HiC-Pro 3.1.0) ----
export HIC_FILE="/path/to/K1_merged.hic"
export MCCOOL="/path/to/K1_merged.mcool"
export SCAFFOLD_LENGTHS="${HI_C_DIR}/data/scaffold.lengths.txt"  # two columns: <scaffold> <length>

# ---- Track visualization inputs ----
export CHROMHMM_SEGMENTS="/path/to/CriGri_PUMC_clean_10_segments.bed"  # chromHMM 10-state
export RNA_BW_CRO="/path/to/Cro.sorted.newref.bigwig"
export RNA_BW_A485="/path/to/10uM.sorted.newref.bigwig"
export A485_RNA_DIR="${RAW_DIR}/A485_CHO_RNAseq/Cleandata"

# ---- onTAD / FAN-C / bedtools ----
export ONTAD_BIN="/path/to/OnTAD"                                 # OnTAD v1.4
export FANC_BIN="/path/to/fanc"                                   # FAN-C v0.9.28
export BEDTOOLS="bedtools"

# ---- 00_hic-pro inputs and helper programs ----
export HIC_PREFIX="CriGri_PUMC"                                       # genome prefix used by 00_hic-pro
export HIC_REF_DIR="/path/to/ref_genome"                          # clean genome, sizes, restriction BED
export HIC_RAW_DIR="${RAW_DIR}/hic_raw"                           # <sample>/<sample>_R{1,2}.fq.gz
export HIC_CLEAN_DIR="${WORK_DIR}/hi_c/clean_fastq"               # trimmed reads
export HICPRO_CONFIG="/path/to/config-CriGri_PUMC_Hae3-bl-hic.txt"
export HICPRO_TRIM_MK="/path/to/BLHiC_filter_trim.mk"
export HICPRO_MAT_MK="/path/to/mat_format.mk"
export HICPRO_ANNO_DIR="/path/to/hicpro_anno_scripts"             # step1_all_chr_size.pl, pureLength_fasta.py
export HICPRO_HOME="/path/to/HiC-Pro_3.1.0"                       # for bin/utils/digest_genome.py

# ==============================================================================
# 07_run_EpiTAD-3D — 3D chromatin reconstruction
# ==============================================================================
export EPITAD3D_HOME="/path/to/EpiTAD-3D"                         # git clone of the EpiTAD-3D repository
export EPITAD3D_MCOOL="${MCCOOL}"                                 # .mcool passed to plot_singleLocus.py
export EPITAD3D_OUT_DIR="${HI_C_DIR}/EpiTAD-3D"                   # output root for the selected loci

# ---- conda environments (names used by the shell scripts) ----
# hifiasm-env          : hifiasm
# nextpolish2-env      : NextPolish2 / meryl / winnowmap / Merqury
# purge_dups-env       : purge_dups / minimap2
# yahs-env             : YaHS / juicer / samtools / BWA-MEM
# gapfill-env          : Flye gap-filling
# busco-env            : BUSCO
# verkko-env           : samtools / pigz (HiFi BAM -> FASTQ)
# tgsgapcloser-env     : minimap2 (ONT BAM construction)
# HiTE                 : de novo transposon prediction, RepeatModeler / RepeatMasker
# repeatmasker-env     : RepeatMasker 4.2.4
# annotation-env       : gene / ncRNA annotation, HISAT2 / samtools
# eviann               : EviAnn
# RNA-seq              : fastp / HISAT2 / samtools / gffread
# CHIP-seq             : fastp / bowtie2 / samtools / MACS3 / deepTools / bedGraphToBigWig
# chromhmm-env         : ChromHMM (Java) / gtfToGenePred / state plots
# hic                  : HiC-Pro deps, valid pairs -> .hic/.mcool, cooler / trackc
# hic-tools            : bedtools / openjdk (Juicer Tools, OnTAD)
# hic-viz-env          : trackC / pyBigWig / cooler track plots
# fanc-env             : FAN-C (KR balance, A/B compartments)
# pastis               : EpiTAD-3D / PASTIS-nb 3D reconstruction
# mummer-env           : MUMmer4 (nucmer)
# genomesyn            : GenomeSyn
# quartet              : quarTeT / CentroMiner
# kmc-env              : KMC
# deeptools-env        : deepTools bamCoverage

# Initialize conda (for `conda activate` in the scripts)
eval "$(conda shell.bash hook)" 2>/dev/null || true
