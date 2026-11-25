# 03_pur — PUR / synteny analysis

| Script | Purpose |
|---|---|
| `run_nucmer.sh` | MUMmer4 nucmer whole-genome alignment (vs CriGri-PICRH-1.0) |
| `pur_stats.py` | PUR statistics (covered / uncovered intervals) |
| `pur_dominant.py` | dominant-chromosome PUR analysis (strict definition) |
| `synteny.py` | synteny block analysis |
| `blast_pur.sh`, `blast_stats.py` | PUR cross-species BLAST conservation (mouse / rat) |
| `pur_gene_classification.py` | PUR-associated gene classification |
| `chr_assignment.awk`, `chr_map.awk` | chromosome assignment / name mapping |

## Environment

| Environment | Purpose |
|---|---|
| `mummer-env` | MUMmer4 (`nucmer`, `mummerplot`) |
| `busco-env` | BLAST+ (`blastn`, `makeblastdb`) for the cross-species BLAST |

Environment files are in `../envs/`.

Software versions: see the "Software versions" section of the top-level `README.md`.
