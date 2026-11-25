# 01_assembly — genome assembly

Pipeline: `HiFi assembly -> polish -> deduplicate -> Hi-C scaffolding -> manual curation -> gap filling -> final polish -> QC`

| Step | Script | Input -> Output | Tool |
|---|---|---|---|
| 0a | `00_prep/convert_hifi_bam2fastq.sh` | HiFi BAM -> `CH6F_G1-{1,2}_hifi.fq.gz` | samtools, pigz |
| 0b | `00_prep/run_fastp_yak.sh` | NGS + Hi-C raw reads -> clean reads + yak k21/k31 | fastp, yak |
| 1 | `01_hifiasm/run.sh` | HiFi -> `CH6F_G1_hifiasm.bp.p_ctg.fasta` | hifiasm 0.25.0 |
| 2 | `02_polish1/run.sh` | p_ctg + HiFi + NGS -> `cho_hifiasm_v3_polished1.fa` | NextPolish2 0.2.2 |
| 3 | `03_purge_dups/run.sh` | polished1 + HiFi -> `purged_clean.fa` | purge_dups |
| 4 | `04_yahs/run.sh` | purged + Hi-C -> `cho_yahs_v3_scaffolds_final.fa` | YaHS 1.2.2 |
| 4b | `04_yahs/make_hic.sh` | scaffolds + BAM -> `.hic` (Juicebox) | juicer_tools |
| 5 | `05_juicebox/curation.sh` | reviewed assembly -> `cho_v3_review4.fa` | juicer post |
| 6a | `06_gapfill/build_bam.sh` | review4 + ONT -> `ont_to_review4.bam` | minimap2 |
| 6b | `06_gapfill/run.sh` | review4 + BAM -> `cho_final_flye_gapfilled.fa` | Flye 2.9.6 |
| 6c | `06_gapfill/fill_gaps_flye.py` | core gap-filling logic | Flye / minimap2 / samtools |
| 7 | `07_polish2/run.sh` | gapfilled + HiFi + NGS -> `cho_review4_gapfilled_polished.fa` | NextPolish2 |
| 8a | `08_qc/run_busco.sh` | genome -> BUSCO (mammalia_odb10) | BUSCO 6.0.0 |
| 8b | `08_qc/run_merqury.sh` | genome + NGS -> QV / k-mer completeness | Merqury 1.4.1 |
| 8c | `08_qc/run_tidk.sh` | genome -> telomere TTAGGG | tidk 0.2.65 |
| 8d | `08_qc/genomescope.R` | KMC3 k-mer histogram -> genome size | GenomeScope 1.0 |
| 8e | centromere prediction | `quarTeT CentroMiner -i <genome> --minperiod 100 --maxperiod 200 --maxgap 50000 --minlength 100000 --identity 0.8 --noplot` | quarTeT 1.3.1 |

## Assembly intermediate files

| Intermediate file | Meaning | Consumed by |
|---|---|---|
| `CH6F_G1_hifiasm.bp.p_ctg.fasta` | hifiasm primary contigs (439, 2.70 Gb) | 02_polish1 |
| `cho_hifiasm_v3_polished1.fa` | first polish | 03_purge_dups |
| `purged_clean.fa` | haplotig-removed (119 contigs) | 04_yahs |
| `cho_yahs_v3_scaffolds_final.fa/.agp` | Hi-C scaffolding (74 scaffolds) | 05_juicebox |
| `cho_v3_review4.fa` | Juicebox manual curation (round 4) | 06_gapfill |
| `cho_final_flye_gapfilled.fa` | Flye gap filling | 07_polish2 |
| `cho_review4_gapfilled_polished.fa` | final polish (review4) | translocation analysis |
| `cho_review5_gapfilled_polished.fa` | **final genome (CriGri_PUMC)** | release |

## Environment

| Environment | Purpose |
|---|---|
| `hifiasm-env` | hifiasm, tidk, gfatools |
| `nextpolish2-env` | NextPolish2, meryl, winnowmap, yak, Merqury, samtools |
| `purge_dups-env` | purge_dups, minimap2 |
| `yahs-env` | YaHS, bwa, minimap2, samtools |
| `gapfill-env` | Flye, minimap2, samtools |
| `tgsgapcloser-env` | TGS-GapCloser, racon, samtools |
| `verkko-env` | Verkko, samtools, fastp, minimap2 |
| `busco-env` | BUSCO |
| `kmc-env` | KMC |
| `quartet` | quarTeT / CentroMiner, BLAST+, TRF |

Environment files are in `../envs/`.

Software versions: see the "Software versions" section of the top-level `README.md`.
