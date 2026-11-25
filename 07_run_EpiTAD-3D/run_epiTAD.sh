#!/usr/bin/env bash
# env: pastis
# REQUIRES (set in ../config.sh):
#   EPITAD3D_HOME     EpiTAD-3D clone (bin/visualization/plot_singleLocus.py)
#   EPITAD3D_MCOOL    .mcool matrix
#   EPITAD3D_OUT_DIR  output root for the selected loci
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && cd ../ && pwd)/config.sh"

mcool_file="${EPITAD3D_MCOOL}"
resolution=10000
plot_script="${EPITAD3D_HOME}/bin/visualization/plot_singleLocus.py"
out_root="${EPITAD3D_OUT_DIR}"

# LOC_00017606 scaffold_7	95937617	95939429
# scaffold_7 96783600-96784800
# 95437617 97284800
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00017606" \
    -r ${resolution} \
    --chr scaffold_7 -s 95437617 -e 97284800 \
    --max_mem 100 --cores 20 \
    --highlight_start 95937617 --highlight_end 95939429 \
    ${mcool_file}

# LOC_00014741 scaffold_5	174507294	174575958
# scaffold_5 174403200-174403400
# 173903200 175075958
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00014741" \
    -r ${resolution} \
    --chr scaffold_5 -s 173903200 -e 175075958 \
    --max_mem 100 --cores 20 \
    --highlight_start 174507294 --highlight_end 174575958 \
    ${mcool_file}


# LOC_00010602 scaffold_3	240425914	240440747
# scaffold_3 240314000-240315800
# 239814000 240940747
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00010602" \
    -r ${resolution} \
    --chr scaffold_3 -s 239814000 -e 240940747 \
    --max_mem 100 --cores 20 \
    --highlight_start 240425914 --highlight_end 240440747 \
    ${mcool_file}

# LOC_00014148 scaffold_5	97354995	97357297
# scaffold_5 97288400-97288600
# 96788400 97857297
# LOC_00014148
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00014148" \
    -r ${resolution} \
    --chr scaffold_5 -s 96788400 -e 97857297 \
    --max_mem 100 --cores 20 \
    --highlight_start 97354995 --highlight_end 97357297 \
    ${mcool_file}


# LOC_00003154
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00003154" \
    -r ${resolution} \
    --chr scaffold_1 -s 452334600 -e 453387424 \
    --max_mem 100 --cores 20 \
    --highlight_start 452878826 --highlight_end 452887424 \
    ${mcool_file}

# LOC_00018041
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00018041" \
    -r ${resolution} \
    --chr scaffold_7 -s 112127000 -e 118723000 \
    --max_mem 100 --cores 20 \
    --highlight_start 117476241 --highlight_end 117499604 \
    ${mcool_file}

# LOC_00018045
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00018045" \
    -r ${resolution} \
    --chr scaffold_7 -s 112130800 -e 118729200 \
    --max_mem 100 --cores 20 \
    --highlight_start 117769246 --highlight_end 117788038 \
    ${mcool_file}

# LOC_00018034
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00018034" \
    -r ${resolution} \
    --chr scaffold_7 -s 114933200 -e 118268600 \
    --max_mem 100 --cores 20 \
    --highlight_start 117082828 --highlight_end 117122241 \
    ${mcool_file}

# LOC_00006965
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00006965" \
    -r ${resolution} \
    --chr scaffold_2 -s 296680800 -e 299411000 \
    --max_mem 100 --cores 20 \
    --highlight_start 298500467 --highlight_end 298856838 \
    ${mcool_file}



# LOC_00010101 scaffold_3 175804101-175807458
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00010101" \
    -r ${resolution} \
    --chr scaffold_3 -s 174680000 -e 176885000 \
    --max_mem 100 --cores 20 \
    --highlight_start 175804101 --highlight_end 175807458 \
    ${mcool_file}

# LOC_00001648 scaffold_1 235045653-235052420
python "${plot_script}" \
    -c CriGri_PUMC -o "${out_root}/LOC_00001648" \
    -r ${resolution} \
    --chr scaffold_1 -s 233730000 -e 237550000 \
    --max_mem 100 --cores 20 \
    --highlight_start 235045653 --highlight_end 235052420 \
    ${mcool_file}
