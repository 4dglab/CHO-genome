# From 1D to 3D Genome of the Chinese Hamster Ovary Cell
Reproducible code for the Chinese hamster ovary (CHO) genome study reported in
*From 1D to 3D Genome of the Chinese Hamster Ovary Cell*.
The repository covers de novo genome assembly, genome annotation, PUR/synteny analysis,
GenomeSyn visualization, chromatin-state segmentation (ChromHMM), Hi-C analysis
(loops, nested TADs, A/B compartments) and 3D chromatin reconstruction.

Final assembly: **CriGri_PUMC** (74 scaffolds, 2.67 Gb total span, ~2.41 Gb haploid genome size).

---

## Contents

- [From 1D to 3D Genome of the Chinese Hamster Ovary Cell](#from-1d-to-3d-genome-of-the-chinese-hamster-ovary-cell)
  - [Contents](#contents)
  - [Data availability](#data-availability)
  - [Requirements](#requirements)
  - [Directory structure](#directory-structure)
  - [Quick start](#quick-start)
  - [Module documentation](#module-documentation)
  - [Conda environments](#conda-environments)
  - [Software versions](#software-versions)
    - [1. Genome assembly](#1-genome-assembly)
    - [2. Assembly QC](#2-assembly-qc)
    - [3. Repeat annotation](#3-repeat-annotation)
    - [4. Gene annotation](#4-gene-annotation)
    - [5. ChIP-seq and RNA-seq preprocessing](#5-chip-seq-and-rna-seq-preprocessing)
    - [6. Chromatin states (ChromHMM)](#6-chromatin-states-chromhmm)
    - [7. Hi-C and 3D genome analysis](#7-hi-c-and-3d-genome-analysis)
    - [8. Statistics and scripting environment](#8-statistics-and-scripting-environment)
    - [Version notes](#version-notes)
  - [Notes and limitations](#notes-and-limitations)
  - [Citation and license](#citation-and-license)
  - [Contact](#contact)

---

## Data availability

| Data | Location |
|---|---|
| ChIP-seq (public) | NCBI SRA PRJNA865478; https://cho-epigenome.boku.ac.at/download.html |
| Control RNA-seq (public) | NCBI SRA PRJNA1123849 |
| Hi-C matrices (`.hic` / `.mcool`) | produced from GSA PRJCA028069 with HiC-Pro 3.1.0; see `06_Hi-C/00_hic-pro/` |
| EpiTAD-3D (3D reconstruction code) | https://github.com/SaotomeMirai/EpiTAD-3D |

> **Hi-C matrices**: the `.hic` / `.mcool` matrices analysed in `06_Hi-C/` were produced from the
> Hi-C raw data (GSA PRJCA028069) with HiC-Pro 3.1.0. Set `HIC_FILE` / `MCCOOL` to these matrices
> before running the Hi-C scripts (see [`06_Hi-C/README.md`](06_Hi-C/README.md)).

> **Raw data filenames**: edit `config.sh` to point each raw-data variable to the file you
> downloaded from GSA. The variable names (`HIFI_BAM_1`, `ONT_UL_1`, `HIC_LIB1_R1`, `NGS_R1`, …)
> are documented in `config.sh`.

---

## Requirements

- Linux x86_64, `bash` ≥ 4, `make`, `awk`, `sort`, `java`
- Conda or mamba (`conda env create -f envs/<name>.yml`)
- CPU/thread count and RAM: set `THREADS` in `config.sh` (default 48). Hi-C and ChromHMM steps use
  large heaps (`-Xms80g -Xmx80g` in `05_chromHMM/01_chromHMM.sh`); adjust to the available RAM.

Programs that are **not** available from conda and must be obtained separately:

| Program | Version | Used by | Configuration |
|---|---|---|---|
| Juicer Tools | 1.22.01 | `06_Hi-C/00_hic-pro/`, `06_Hi-C/01_loop_calling/` | `JUICER_TOOLS` in `config.sh` |
| OnTAD | v1.4 | `06_Hi-C/02_tad/` | `ONTAD_BIN` |
| HiC-Pro | 3.1.0 | `06_Hi-C/00_hic-pro/` | source install; `HiC-Pro` on `PATH` |
| ChromHMM | 1.26 | `05_chromHMM/` | `chromHMM_file` in `01_chromHMM.sh` |
| GenomeSyn | 1.2.7 | `04_pur_visualization/` | `GENOMESYN` in `config.sh` |
| gnuplot | — | `03_pur/` (`mummerplot`) | `GNUPLOT_BIN`, `GNUPLOT_PS_DIR` |
| EpiTAD-3D | v1.0 | `07_run_EpiTAD-3D/` | `git clone https://github.com/SaotomeMirai/EpiTAD-3D` |

---

## Directory structure

```
CHO_publish/
├── README.md                 # this file (global setup, environments, software versions)
├── config.sh                 # global configuration (edit this first)
├── envs/                     # 25 conda environment files
├── 01_assembly/              # genome assembly
│   └── README.md             # module documentation
├── 02_annotation/            # genome annotation + ChIP-seq / RNA-seq alignment rate
│   └── README.md
├── 03_pur/                   # PUR / synteny analysis
│   └── README.md
├── 04_pur_visualization/     # GenomeSyn visualization
│   └── README.md
├── 05_chromHMM/              # chromatin-state segmentation
│   └── README.md
├── 06_Hi-C/                  # Hi-C processing and 3D genome analysis
│   └── README.md
└── 07_run_EpiTAD-3D/         # 3D chromatin reconstruction (EpiTAD-3D)
    └── README.md
```

Each analysis module has its own `README.md`. Script-by-script documentation, inputs,
outputs and per-module environments are documented there.

---

## Quick start

1. **Create the environments.** One file per environment:
   `conda env create -f envs/hifiasm-env.yml` (repeat for the environments you need).
   The name inside each file is the name used by the scripts.
2. **Edit `config.sh`.** Set `RAW_DIR`, `WORK_DIR`, `THREADS`, `CONDA_BASE`, `JUICER_TOOLS`,
   `BUSCO_DOWNLOAD`, `GENOMESYN`, the gnuplot paths, the reference-genome paths and the
   module blocks (`02_annotation`, `05_chromHMM`, `06_Hi-C`, `07_run_EpiTAD-3D`).
3. **Source `config.sh` in the shell.** The `.sh` scripts source it themselves. The Snakemake
   workflows and the Python scripts read the same settings from the environment, so run
   `source config.sh` before them.
4. **Run the modules in order**: `01_assembly` → `02_annotation` → `03_pur` →
   `04_pur_visualization` → `05_chromHMM` → `06_Hi-C` → `07_run_EpiTAD-3D`.

Scripts marked `set -euo pipefail` read configuration through
`source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"`.

---

## Module documentation

Script-by-script documentation is in the `README.md` of the corresponding module directory.

| Module | Documentation | Scope |
|---|---|---|
| 01_assembly | [`01_assembly/README.md`](01_assembly/README.md) | HiFi assembly, polishing, Hi-C scaffolding, gap filling, assembly QC |
| 02_annotation | [`02_annotation/README.md`](02_annotation/README.md) | repeat, gene and ncRNA annotation, ChIP-seq / RNA-seq alignment rate |
| 03_pur | [`03_pur/README.md`](03_pur/README.md) | PUR, synteny and cross-species BLAST |
| 04_pur_visualization | [`04_pur_visualization/README.md`](04_pur_visualization/README.md) | GenomeSyn synteny, TE and gene-density plots |
| 05_chromHMM | [`05_chromHMM/README.md`](05_chromHMM/README.md) | chromatin-state segmentation (10 states) |
| 06_Hi-C | [`06_Hi-C/README.md`](06_Hi-C/README.md) | HiC-Pro processing, loops, nested TADs, A/B compartments, tracks |
| 07_run_EpiTAD-3D | [`07_run_EpiTAD-3D/README.md`](07_run_EpiTAD-3D/README.md) | 3D chromatin reconstruction |

The order of the modules follows the analysis order:
`01_assembly` → `02_annotation` → `03_pur` → `04_pur_visualization` → `05_chromHMM` →
`06_Hi-C` → `07_run_EpiTAD-3D`.

---

## Conda environments

Each file in `envs/` creates one environment. The name inside the file is the name used by the scripts.
Versions are intentionally not pinned; see [Software versions](#software-versions) for the versions
used in the study.

| Environment | Purpose |
|---|---|
| `hifiasm-env` | hifiasm, tidk, gfatools |
| `nextpolish2-env` | NextPolish2, meryl, winnowmap, yak, Merqury, samtools |
| `purge_dups-env` | purge_dups, minimap2 |
| `yahs-env` | YaHS, bwa, minimap2, samtools |
| `gapfill-env` | Flye gap filling, minimap2, samtools, NextDenovo |
| `tgsgapcloser-env` | TGS-GapCloser, racon, samtools |
| `verkko-env` | Verkko, samtools, fastp, minimap2 |
| `busco-env` | BUSCO |
| `kmc-env` | KMC (k-mer histogram for GenomeScope) |
| `quartet` | quarTeT / CentroMiner, BLAST+, TRF |
| `HiTE` | HiTE, RepeatModeler, RepeatMasker 4.1.1, HISAT2, minimap2 |
| `repeatmasker-env` | RepeatMasker 4.2.4 |
| `annotation-env` | HISAT2, StringTie, Augustus, SNAP, Infernal, tRNAscan-SE, AGAT, eggNOG-mapper, MAKER |
| `eviann` | EviAnn, SNAP |
| `RNA-seq` | fastp, HISAT2, samtools, gffread |
| `CHIP-seq` | fastp, FastQC, bowtie2, samtools, MACS3, deepTools, bedGraphToBigWig, Snakemake |
| `chromhmm-env` | OpenJDK, gtfToGenePred, pandas, matplotlib, seaborn |
| `hic` | HiC-Pro dependencies, cooler, hic2cool, cooltools, bx-python, iced, samtools, bowtie2, R, OpenJDK, Snakemake |
| `hic-tools` | bedtools, OpenJDK (Juicer Tools, OnTAD) |
| `hic-viz-env` | trackC, pyBigWig, cooler, numpy, pandas, matplotlib, scipy |
| `fanc-env` | FAN-C |
| `pastis` | EpiTAD-3D / PASTIS-nb |
| `deeptools-env` | deepTools `bamCoverage` |
| `mummer-env` | MUMmer4 (`nucmer`, `mummerplot`) |
| `genomesyn` | GenomeSyn 1.2.7 (bioconda `genomesyn2`, provides the Perl runtime) |

---

## Software versions

Versions were read from (i) the conda environments used for the analysis, (ii) the pins in
`envs/*.yml`, and (iii) `software_versions.md` of the earlier code release for programs that are not
installed through conda. See [Version notes](#version-notes) for deviations from the earlier list.

### 1. Genome assembly

| Software | Version | Purpose | Environment |
|---|---|---|---|
| hifiasm | 0.25.0 (r726) | HiFi + ONT primary-contig assembly | `hifiasm-env` |
| NextPolish2 | 0.2.2 | two-round polishing with NGS short reads (Polish_1, Polish_2) | `nextpolish2-env` |
| purge_dups | — | haplotypic redundancy removal (439 -> 119 contigs) | `purge_dups-env` |
| YaHS | 1.2.2 | Hi-C chromosome scaffolding | `yahs-env` |
| Juicebox | — | Hi-C scaffolding manual curation (reviews 1-5) | — |
| Flye | 2.9.6 (b1802) | ONT gap filling | `gapfill-env` |
| MUMmer4 (nucmer) | 4.0.1 | PUR / synteny all-to-all alignment | `mummer-env` |
| minimap2 | 2.31 (r1302) | read alignment (HiFi / ONT / short reads) | `nextpolish2-env` |
| samtools | 1.24 | BAM processing | `nextpolish2-env` |
| bedtools | 2.31.1 | interval operations | `nextpolish2-env` |
| meryl | 1.4.2 | k-mer counting (Merqury, NextPolish2) | `nextpolish2-env` |
| fastp | 0.24.0 | Hi-C / NGS read QC and trimming | `verkko-env`, `nextpolish2-env` |
| BWA-MEM | 0.7.19 | Hi-C read alignment (`-5SP`) | `yahs-env`, `verkko-env` |
| yak | — | k-mer databases for NextPolish2 | `nextpolish2-env` |

### 2. Assembly QC

| Software | Version | Purpose | Environment |
|---|---|---|---|
| BUSCO | 6.0.0 | assembly / annotation completeness (mammalia_odb10) | `busco-env` |
| Merqury | 1.4.1 | QV / completeness / spectra (QV 62.85, completeness 98.89%) | `nextpolish2-env` |
| GenomeScope | 1.0 (k=21) | k-mer genome size estimation (2.41 Gb) | — |
| tidk | 0.2.65 | telomere repeat (TTAGGG) identification | `hifiasm-env` |
| quarTeT (CentroMiner) | 1.3.1 | centromere prediction | `quartet` |
| KMC | 3.2.4 | k-mer histogram for GenomeScope | `kmc-env` |
| gap statistics | — | 43 gaps / 3,915 bp (custom scripts) | — |

### 3. Repeat annotation

| Software | Version | Purpose | Environment |
|---|---|---|---|
| HiTE | 3.3.3 | de novo transposon prediction | `HiTE` |
| Dfam | 4.0 | homology-based repeat database | — |
| RepeatMasker | 4.1.1 / 4.2.4 | repeat masking | `HiTE` / `repeatmasker-env` |
| RepeatModeler | 2.0.1 | de novo repeat modelling | `HiTE` |
| TRF | 4.10.0rc2 | tandem repeats | `HiTE` |

### 4. Gene annotation

| Software | Version | Purpose | Environment |
|---|---|---|---|
| HISAT2 | 2.2.3 | RNA-seq alignment | `annotation-env` / `eviann` |
| StringTie | 2.2.3 (eviann: 3.0.0) | transcript assembly | `annotation-env` |
| TransDecoder | — | ORF translation | — |
| Augustus | 3.5.0 | ab initio gene prediction | `annotation-env` |
| SNAP | — | ab initio gene prediction | `annotation-env` |
| EviAnn | 2.0.6 | evidence-driven annotation | `eviann` |
| EggNOG-mapper | 2.1.12 | functional annotation (GO / KEGG / COG) | `annotation-env` |
| BLAST+ (blastn, makeblastdb) | 2.14.1 (quartet: 2.17.0) | novel-gene comparison, PUR cross-species BLAST | `annotation-env`, `busco-env` |
| AGAT | — | GFF processing / filtering | `annotation-env` |
| tRNAscan-SE | 2.0.13 | tRNA prediction | `annotation-env` |
| Infernal / Rfam | 1.1.5 | ncRNA (rRNA / snRNA / miRNA) annotation | `annotation-env` |

### 5. ChIP-seq and RNA-seq preprocessing

| Software | Version | Purpose | Environment |
|---|---|---|---|
| fastp | 0.24.0 | read QC and adapter trimming | `CHIP-seq`, `RNA-seq` |
| FastQC | 0.12.1 | read QC report | `CHIP-seq` |
| Bowtie2 | 2.5.4 | ChIP-seq alignment | `CHIP-seq` |
| HISAT2 | 2.2.1 | RNA-seq alignment and index building | `RNA-seq` |
| samtools | 1.21 | BAM conversion, sorting, indexing | `CHIP-seq`, `RNA-seq` |
| MACS3 | 3.0.2 | peak calling (broad and narrow), FE bedgraph | `CHIP-seq` |
| UCSC bedGraphToBigWig | 482 | bedGraph to bigWig | `CHIP-seq` |
| deepTools (bamCoverage) | 3.5.6 | bigWig generation (RPKM) | `deeptools-env` |
| gffread | — | GFF to GTF conversion for HISAT2 index | `RNA-seq` |
| Snakemake | 7.32.4 / 8.27.1 / 9.13.3 | workflow engine for the `.smk` files | `hic` / `CHIP-seq` / `RNA-seq` |

### 6. Chromatin states (ChromHMM)

| Software | Version | Purpose | Environment |
|---|---|---|---|
| ChromHMM | 1.26 | BinarizeBam, LearnModel (10 states), OverlapEnrichment | `chromhmm-env` (jar) |
| UCSC gtfToGenePred | — | gene table for `ConvertGeneTable` | `chromhmm-env` |
| OpenJDK | 23.0.1 / 23.0.2 | Java runtime for ChromHMM and Juicer Tools | `chromhmm-env`, `CHIP-seq` / `hic` |
| pandas | — | state tables in the notebooks | `chromhmm-env` |
| matplotlib, seaborn | — | state plots and heatmaps | `chromhmm-env` |

### 7. Hi-C and 3D genome analysis

| Software | Version | Purpose | Environment |
|---|---|---|---|
| HiC-Pro | 3.1.0 | BL-Hi-C read processing and filtering | source install (`hic` for its Python deps) |
| Bowtie2 | 2.5.4 | Hi-C read alignment inside HiC-Pro | `hic` |
| samtools | 1.21 | valid-pair and BAM handling | `hic` |
| hic2cool | 1.0.1 | `.hic` to `.mcool` conversion | `hic` |
| cooler | 0.9.3 / 0.10.3 | matrix handling (`.mcool`, `.cool`) | `hic` / `pastis` |
| cooltools | 0.7.1 | matrix analysis | `hic` |
| Juicer Tools (HiCCUPS, pre) | 1.22.01 | loop calling (`-k KR -m 500`) and `.hic` creation | `hic-tools` |
| OnTAD | v1.4 | nested TAD identification (level 1-5) | `hic-tools` |
| FAN-C | 0.9.28 | KR balancing and A/B compartments | `fanc-env` |
| bedtools | 2.31.1 | Hi-C interval operations | `hic-tools` |
| trackC | 0.0.18 | Hi-C / locus track visualization | `hic-viz-env` |
| pyBigWig | 0.3.24 | bigWig track reading | `hic-viz-env`, `hic` |
| OpenJDK | 23.0.2 | Java runtime for Juicer Tools and OnTAD | `hic`, `hic-tools` |
| EpiTAD-3D (Pastis-NB) | v1.0 | 3D structure reconstruction (MDS + NB) | `pastis` |
| PyMOL | 3.1.6.1 | 3D structure visualization | — |
| Cytoscape | 3.10.3 | gene-chromatin interaction network | — |
| HOMER | v5.1 | de novo motif discovery | — |
| JASPAR (CORE vertebrates) | — | TF motif database | — |

### 8. Statistics and scripting environment

| Software | Version | Purpose | Environment |
|---|---|---|---|
| R | 4.4.2 / 4.4.1 / 4.2.3 | statistics and plotting | `hic`, `RNA-seq` / `CHIP-seq` / `HiTE` |
| Python | 3.9.21 / 3.9.23 / 3.11.11 | pipeline scripts | `hic` / `pastis` / `CHIP-seq`, `RNA-seq` |
| GenomeSyn | 1.2.7 | synteny / repeat / gene-density plots | `genomesyn` |
| gnuplot | — | `mummerplot` output for nucmer | — |

### Version notes

Several numbers differ from `software_versions.md` of the earlier release because they were re-read
from the environments used here. The tables above follow the measured values:

| Software | Earlier list | This repository |
|---|---|---|
| fastp | 1.3.3 | 0.24.0 |
| HISAT2 | 2.2.3 | 2.2.1 (RNA-seq alignment); 2.2.3 is the version cited for annotation |
| OpenJDK | 17 | 23.0.1 / 23.0.2 |
| trackC | 0.0.19 | 0.0.18 |
| pyBigWig | 0.3.22 | 0.3.24 |
| cooler | 0.10.2 | 0.9.3 (`hic`), 0.10.3 (`pastis`) |
| R | 4.5.2 | 4.4.2 / 4.4.1 |
| Python | 3.12.7 / 3.14 | 3.9.21 / 3.9.23 / 3.11.11 |

Reconcile these with the manuscript Methods before submission.

---

## Notes and limitations

- **Unpinned environments.** `envs/*.yml` do not pin most versions. For exact reproduction, pin the
  versions listed above.
- **Configuration.** Every path is set in `config.sh`. The `.sh` scripts source it themselves.
  The Snakemake files (`.smk`) and the plotting scripts read the same values from the environment,
  so run `source config.sh` first. Each script lists its inputs in a `# REQUIRES` header.
- **One file to edit by hand.** `06_Hi-C/00_hic-pro/config-CriGri_PUMC_Hae3-bl-hic.txt` is a HiC-Pro
  configuration file and cannot read `config.sh`. Edit its three `/path/to` entries.
- **External helper programs.** `06_Hi-C/00_hic-pro/` needs `step1_all_chr_size.pl`,
  `pureLength_fasta.py`, `digest_genome.py` and the trimming/matrix makefiles. Point
  `HICPRO_ANNO_DIR`, `HICPRO_HOME`, `HICPRO_TRIM_MK` and `HICPRO_MAT_MK` to your copies.
- **Notebooks** (`05_chromHMM/*.ipynb`) read `CHROMHMM_PROJECT_DIR` and `CHROMHMM_GENOME` from the
  environment (defaults: `chromHMM_project`, `CriGri_PUMC_clean`).

---

## Citation and license

This repository is released under the MIT License. See [`LICENSE`](LICENSE).

## Contact

If you have any questions, please open an issue on GitHub or contact:  
**gy1.mirai2003@gmail.com**

---

© 2026 **4DGenome Lab, PUMC**. All rights reserved.
