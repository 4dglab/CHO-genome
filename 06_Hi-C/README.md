# 06_Hi-C — Hi-C / 3D genome

## Settings in `config.sh`

The Hi-C settings are in `../config.sh`. Each script repeats the list in its `# REQUIRES` header:

| Variable | Content |
|---|---|
| `HI_C_DIR` | output root of the Hi-C module (e.g. `${WORK_DIR}/hi_c`) |
| `HIC_FILE` | merged `.hic` matrix from HiC-Pro and `juicer pre` |
| `MCCOOL` | merged `.mcool` matrix from `hic2cool` |
| `GENOME_FA` | reference genome FASTA with `.fai` |
| `GENOME_GFF` | gene annotation GFF3 |
| `GENOME_BED12` | gene annotation in BED12 format |
| `CHROMHMM_SEGMENTS` | chromHMM 10-state segments BED |
| `HISTONE_BW_DIR` | `<mark>/<mark>_FE.bw` directory (from `../05_chromHMM/peak_signal_*.smk`) |
| `RNA_BW_CRO`, `RNA_BW_A485` | RNA-seq bigwigs (control, treated) |
| `A485_RNA_DIR` | clean A485 RNA-seq reads `<sample>/<sample>_R{1,2}.fq.gz` |
| `SCAFFOLD_LENGTHS` | two-column TSV `<scaffold> <length>` (`cut -f1,2 genome.fa.fai`) |
| `ONTAD_BIN` | OnTAD v1.4 executable |
| `FANC_BIN` | `fanc` executable |
| `BEDTOOLS` | `bedtools` executable |
| `HIC_PREFIX` | genome prefix used by `00_hic-pro` |
| `HIC_REF_DIR` | clean genome, chrom sizes and restriction BED for `00_hic-pro` |
| `HIC_RAW_DIR` | raw Hi-C reads `<sample>/<sample>_R{1,2}.fq.gz` |
| `HIC_CLEAN_DIR` | trimmed Hi-C reads (written by the workflow) |
| `HICPRO_CONFIG` | HiC-Pro configuration file |
| `HICPRO_TRIM_MK` | linker-trimming makefile |
| `HICPRO_MAT_MK` | matrix-conversion makefile (`juicer pre`, `hic2cool`) |
| `HICPRO_ANNO_DIR` | `step1_all_chr_size.pl` and `pureLength_fasta.py` |
| `HICPRO_HOME` | HiC-Pro install (`bin/utils/digest_genome.py`) |

`JUICER_TOOLS` and `THREADS` come from `../config.sh` as well.

## 00_hic-pro — BL-Hi-C processing

| Step | Script | Purpose | Tool |
|---|---|---|---|
| 1 | `00_hic-pro/01-hic-pro_annofiles.sh` | clean genome (scaffold_1-11), chrom sizes, bowtie2 index, HaeIII restriction fragments | samtools, bowtie2-build, HiC-Pro `digest_genome.py`, custom Perl |
| 2 | `00_hic-pro/02-run_hic-pro.smk` | linker trimming -> HiC-Pro 3.1.0 -> `allValidPairs` -> `.hic` -> `.mcool` | HiC-Pro, juicer_tools, hic2cool |

### Reference preparation

`*_clean` means `scaffold_1` to `scaffold_11` only. The other 63 scaffolds are
removed before Hi-C and chromatin-state analysis.

`00_hic-pro/01-hic-pro_annofiles.sh` builds the clean set from the full assembly:

1. `samtools faidx` indexes the full FASTA.
2. A loop extracts `scaffold_1` to `scaffold_11` into `<prefix>_clean.fa`.
3. `step1_all_chr_size.pl` writes `<prefix>_clean.sizes`.
4. `bowtie2-build` builds the index.
5. `digest_genome.py` writes the HaeIII fragment table.

The clean GTF is the full GTF restricted to the same 11 scaffolds.
The original command was not recorded. This filter reproduces it:

```bash
awk -F'\t' 'NR<=3 || $1 ~ /^scaffold_([1-9]|1[01])$/' CriGri_PUMC.gtf > CriGri_PUMC_clean.gtf
```

The step is optional. Point `GENOME_FA`, `GENOME_SIZES` and `GENOME_GTF` at the full
assembly if the indices already match it.

The HiC-Pro configuration used in the study is `00_hic-pro/config-CriGri_PUMC_Hae3-bl-hic.txt`
(BL-Hi-C, HaeIII, ligation site `CGCGATATCTTATCTGAC`).
Trimming and matrix conversion use the `BLHiC_filter_trim.mk` and `mat_format.mk` makefiles(in-house written, equivalent software may be used instead.)

## 01_loop_calling

| Step | Script | Input -> Output | Tool |
|---|---|---|---|
| 1 | `01_loop_calling/run_hiccups.sh` | `K1_merged.hic` -> `loop/K1_merged.hiccups/merged_loops.bedpe` | Juicer Tools 1.22.01 (HiCCUPS, `-k KR -m 500`) |

## 02_tad — nested TADs (10 kb)

| Step | Script | Input -> Output | Tool |
|---|---|---|---|
| 1 | `02_tad/01_run_ontad.sh` | `.hic` + `SCAFFOLD_LENGTHS` -> `TAD/10kb.onTAD/onTAD.K1.<scaffold>.{tad,bed}` | OnTAD v1.4 |
| 2 | `02_tad/02_level_label_bed.sh` | OnTAD `.bed` -> level1 to level5 labels | awk |
| 3 | `02_tad/03_encode_tad.py` | `.tad` -> `TAD/combined/onTAD.10kb.encoded.tsv` | Python |
| 4 | `02_tad/04_make_annotation.sh` | encoded TSV -> `TAD/combined/genome.10kb.onTAD.annotation.bed` (gap / Singleton / level1-5) | awk, bedtools |
| 5 | `02_tad/05_category_count.sh` | annotation -> `TAD/combined/onTAD.category_count.10kb.tsv` | awk |

## 03_compartments — A/B (1 mb, 100 kb)

| Step | Script | Input -> Output | Tool |
|---|---|---|---|
| 1 | `03_compartments/01_balance.sh` | `.mcool` slices -> `K1_{1mb,100kb}_balanced.hic` | FAN-C 0.9.28 (`fanc hic -m KR`) |
| 2 | `03_compartments/02_fanc_compartments.sh` | balanced `.hic` + genome -> `K1.{1mb,100kb}_balanced.{ab,bed,ev.txt}` | FAN-C (`fanc compartments`) |
| 3 | `03_compartments/03_tad_ab_intersect.sh` | onTAD 10 kb x A/B bed -> `onTAD.10kb.intersect.AB<RES>.{bed,stats.tsv,level_distribution.tsv}` | bedtools, awk |

## 04_visualization

| Script | Purpose | Tool |
|---|---|---|
| `04_visualization/run_rnaseq_bigwig.sh` | align A485 Cro / 10uM RNA-seq, build RPKM bigwigs | HISAT2, samtools, deepTools `bamCoverage` |
| `04_visualization/make_fig6_tracks.py` | locus tracks: loop links, 6 histone marks, RNA-seq (Cro / A485), gene model, Hi-C triangle, onTAD markline | trackC, pyBigWig |
| `04_visualization/show_region_trackc.py` | region tracks: chromHMM states, histone marks, RNA-seq, gene models | trackC, pyBigWig |

## Environments

| Step | Environment (`../envs/`) | Notes |
|---|---|---|
| HiC-Pro and matrices | `hic` | HiC-Pro 3.1.0 also needs bowtie2, samtools, R, python; install it from source |
| Loop calling | `hic-tools` | Java runtime for the Juicer Tools jar |
| Nested TAD | `hic-tools` | bedtools, Java for OnTAD |
| A/B compartments | `fanc-env` | FAN-C |
| Visualization | `hic-viz-env` | trackC, pyBigWig, cooler, matplotlib |
| RNA-seq bigwig | `annotation-env` + `deeptools-env` | HISAT2 / samtools then bamCoverage |

Software versions: see the "Software versions" section of the top-level `README.md`.
