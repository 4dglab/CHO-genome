#!/usr/bin/env bash
# PURPOSE: Intersect the 10 kb onTAD annotation with A/B compartments (1 mb and 100 kb)
#          and summarize onTAD-level distributions within A/B.
# INPUT:   ${HI_C_DIR}/TAD/combined/genome.10kb.onTAD.annotation.bed
#          ${HI_C_DIR}/AB_compartments/K1.{1mb,100kb}_balanced.bed
# OUTPUT (per resolution, into ${HI_C_DIR}/TAD/combined/):
#          onTAD.10kb.intersect.AB<RES>.bed                 (bedtools -f 0.7 -wo, dedup)
#          onTAD.10kb.intersect.AB<RES>.stats.tsv           (A/B x level counts)
#          onTAD.10kb.intersect.AB<RES>.level_distribution.tsv (per-level A/B %)
# TOOL:    bedtools / awk
# env: hic-tools
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   BEDTOOLS           bedtools executable (hic-tools env)
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../../ && pwd)/config.sh"

ONTAD="${HI_C_DIR}/TAD/combined/genome.10kb.onTAD.annotation.bed"
ABDIR="${HI_C_DIR}/AB_compartments"
OUTDIR="${HI_C_DIR}/TAD/combined"
mkdir -p "${OUTDIR}"

for RES in 1mb 100kb; do
    AB="${ABDIR}/K1.${RES}_balanced.bed"
    OUT="${OUTDIR}/onTAD.10kb.intersect.AB${RES}.bed"
    STATS="${OUTDIR}/onTAD.10kb.intersect.AB${RES}.stats.tsv"
    DIST="${OUTDIR}/onTAD.10kb.intersect.AB${RES}.level_distribution.tsv"

    echo "===== $(date) onTAD(10kb) x AB(${RES}) ====="
    "${BEDTOOLS}" intersect -a "${ONTAD}" -b "${AB}" -f 0.7 -wo > "${OUT}"
    sort -u "${OUT}" -o "${OUT}"

    awk -F'\t' '{
        group=$8; level=$4
        count[group,level]++; total[group]++
    }
    END {
        split("gap Singleton level1 level2 level3 level4 level5", levels, " ")
        printf "Group\tTotal\t"
        for (i in levels) printf "%s\t", levels[i]
        print ""
        for (group in total) {
            printf "%s\t%d\t", group, total[group]
            for (i in levels) { level=levels[i]; v=((group SUBSEP level) in count)?count[group,level]:0; printf "%d\t", v }
            print ""
        }
    }' "${OUT}" > "${STATS}"
    cat "${STATS}"

    awk -F'\t' '
      { c[$4]++; ab[$4,$8]++ }
    END {
      n=split("gap Singleton level1 level2 level3 level4 level5", levels, " ")
      printf "level\tA\tB\ttotal\tA_pct\tB_pct\n"
      for (i=1;i<=n;i++){ l=levels[i]; a=ab[l,"A"]; b=ab[l,"B"]; t=a+b
        ap=(t==0)?0:a/t*100; bp=(t==0)?0:b/t*100
        printf "%s\t%d\t%d\t%d\t%.2f\t%.2f\n", l,a,b,t,ap,bp }
    }' "${OUT}" > "${DIST}"
    cat "${DIST}"
    echo "DONE ${RES}: ${OUT}"
done
echo "ALL_TAD_AB_DONE"
