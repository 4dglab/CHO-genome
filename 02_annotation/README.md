# 02_annotation — genome annotation

Pipeline: `repeat annotation (HiTE + Dfam) -> gene structure (Augustus + EviAnn + AGAT) -> filtering -> ncRNA -> QC`

| Step | Script | Tool |
|---|---|---|
| A | `stepA_repeat/run_hite.sh`, `run_repeatmasker.sh`, `run_merge_softmask.sh` | HiTE 3.3.3 / RepeatMasker 4.2.4 |
| B | `stepB_gene/run_stepB_all.sh` (B0 prep -> B1 Augustus -> B2 EviAnn -> B3 AGAT) | HISAT2 / StringTie / Augustus / EviAnn / AGAT |
| B4 | `stepB4_clean/run_stepB4_clean.sh` (redundancy removal, pseudogene filtering, evidence tiers L1-L4) | custom |
| C | `stepC_ncrna/run_stepC_ncrna.sh` | Infernal cmscan / tRNAscan-SE 2.0 |
| D | `stepD_qc/run_stepD_busco_eggnog.sh` | BUSCO 6.0.0 / eggNOG-mapper 2.1.12 |
| — | `scripts/` (`augustus_build_training.py`, `augustus_chunk_predict.py`, `evidence_classify.py`, `filtering_report.py`, `gtf2gff3.py`, `merge_repeat_gff3.py`, `overlap_analysis.py`, `pseudogene_scan.py`, `softmask_genome.py`) | custom Python helpers |

## Alignment-rate comparison (`alignment_rate/`)

| Script | Purpose | Environment |
|---|---|---|
| `ChIP-seq/01-CHIP_mapping_rate.smk` | ChIP-seq read QC and mapping rate against each genome (cg20, cg1001, cho_ont) | `CHIP-seq` |
| `RNA-seq/02-RNA-seq_mapping_rate.smk` | RNA-seq read QC and mapping rate against each genome | `RNA-seq` |
| `RNA-seq/03-hisat2-build.sh` | HISAT2 index with splice sites and exons from GTF | `RNA-seq` |

These scripts read `ANNOT_REF_DIR`, `ONT_REF_DIR`, `BT2_INDEX_DIR`, `HISAT2_INDEX_DIR`,
`CHIP_ALIGN_DIR` and `RNA_ALIGN_DIR` from `../config.sh`.
Run `source config.sh` before starting the Snakemake workflows.

## Environment

| Environment | Purpose |
|---|---|
| `HiTE` | HiTE, RepeatModeler, RepeatMasker 4.1.1, HISAT2 |
| `repeatmasker-env` | RepeatMasker 4.2.4 |
| `annotation-env` | HISAT2, StringTie, Augustus, SNAP, Infernal, tRNAscan-SE, AGAT, eggNOG-mapper |
| `eviann` | EviAnn, SNAP |
| `busco-env` | BUSCO |
| `CHIP-seq` | fastp, bowtie2, samtools (alignment-rate comparison) |
| `RNA-seq` | fastp, HISAT2, samtools, gffread (alignment-rate comparison) |

Environment files are in `../envs/`.

Software versions: see the "Software versions" section of the top-level `README.md`.
