#!/usr/bin/env bash
# PURPOSE: RepeatMasker (Dfam) homology-based repeat annotation
# INPUT:   genome fasta
# OUTPUT:  genome.fa.out
# TOOL:    RepeatMasker 4.2.4 (-species Cricetulus griseus)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"


# Step 1.5: RepeatMasker with Dfam 4.0 (formal annotation)
set -eo pipefail
conda activate repeatmasker-env

BASE=${WORK_DIR}/annotation
cd ${BASE}/04_repeatmasker

echo "[$(date)] RepeatMasker (Dfam 4.0, Cricetulus griseus) start"
RepeatMasker -species "Cricetulus griseus" -pa 32 -gff -dir . ${BASE}/01_genome/genome.fa 2>&1 | tee ${BASE}/09_logs/repeatmasker_dfam.log
echo "[$(date)] RepeatMasker done"

echo "[$(date)] Summary"
if [ -f genome.fa.out ]; then
    echo "Masking summary (genome.fa.out):"
    tail -5 genome.fa.out
fi
