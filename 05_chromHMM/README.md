# 05_chromHMM — chromatin states

Pipeline: `ChIP-seq BAM -> bigwig -> BinarizeBam -> LearnModel (10 states) -> OverlapEnrichment`.

| Step | Script | Input -> Output | Tool |
|---|---|---|---|
| 1 | `00_oribam2bam.smk` | original ChIP-seq BAM -> fastq -> bowtie2 -> `sorted_bam/*_sorted.bam(.csi)` | samtools, bowtie2 |
| 2 | `peak_signal_narrow.smk` | `H3K4me3`, `H3K27ac` -> MACS3 narrow peaks -> FE bedgraph -> `bw/<mark>/<mark>_FE.bw` | MACS3, bedGraphToBigWig |
| 3 | `peak_signal_broad.smk` | `H3K4me1`, `H3K9me3`, `H3K27me3`, `H3K36me3` -> MACS3 broad peaks -> `bw/<mark>/<mark>_FE.bw` | MACS3, bedGraphToBigWig |
| 4 | `01_chromHMM.sh` | sorted BAMs + `mark_table.tsv` -> `chromHMM_out/model/<genome>_10_segments.bed` | ChromHMM 1.26 |
| 5 | `02_visual.ipynb` | state tracks over selected loci | pandas, matplotlib |
| 6 | `states_heatmap.ipynb` | state heatmaps | pandas, seaborn |

Notes:

- Steps 1 to 3 run under the `CHIP-seq` environment; steps 4 to 6 run under `chromhmm-env`.
- Paths come from `../config.sh`: `GENOME_FA`, `GENOME_SIZES`, `GENOME_GTF`, `CHIP_BAM_DIR`,
  `CHIP_ORI_BAM_DIR`, `HISTONE_BW_DIR`, `MACS3_OUT_DIR`, `CHROMHMM_PROJECT_DIR`, `CHROMHMM_HOME`,
  `CHROMHMM_GENOME`, `GTF_TO_GENEPRED` and `CHROMHMM_JAVA_OPTS`.
  The `.smk` files read these from the environment, so run `source config.sh` first.
- `GENOME_GTF` is the clean annotation (`scaffold_1` to `scaffold_11`).
  See [06_Hi-C/README.md](../06_Hi-C/README.md#reference-preparation) for how it is built.
- `mark_table.tsv` lists the BAM file and the mark name. Paths inside that table must be **relative**,
  because ChromHMM concatenates paths internally.
- ChromHMM also needs a chrom-sizes file in `CHROMHMM/CHROMSIZES/<genome>.txt` and the gene table
  converted with `ConvertGeneTable` (the script calls `gtfToGenePred` first).

## Environment

| Environment | Purpose |
|---|---|
| `CHIP-seq` | bowtie2, samtools, MACS3, bedGraphToBigWig (steps 1-3) |
| `chromhmm-env` | ChromHMM (Java), gtfToGenePred, pandas, matplotlib, seaborn (steps 4-6) |

Environment files are in `../envs/`.

Software versions: see the "Software versions" section of the top-level `README.md`.
