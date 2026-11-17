#!/bin/bash
#SBATCH --job-name=chip_signal_dars
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=bbdefault
#
# Quantify ChIP-seq signal over the regions that gained or lost accessibility in
# patients, and draw profile + heatmap plots (Figure 10 of the report).
#
# Rather than calling ChIP peaks first, the patient DNase regions are the
# reference and the ChIP signal is measured over them (+/- 2 kb around the
# centre). This was our supervisor's suggestion: the question is what binds at
# the regions that change in patients, so those regions define where to look.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2022a
module load deepTools/3.5.2-foss-2022a

mkdir -p "$CHIP_PLOTS"
tracks=(HA-RUNX1-ETO_0Dox HA-RUNX1-ETO_5Dox RUNX1_0Dox RUNX1_5Dox H3K27ac_0Dox H3K27ac_5Dox)

bigwigs=(); labels=()
for t in "${tracks[@]}"; do
    f=$(ls "$CHIP_BW"/*"_${t}.bw" | head -n 1)
    bigwigs+=("$f"); labels+=("${t#HA-}")
done

for set in gained lost; do
    computeMatrix reference-point --referencePoint center \
        -S "${bigwigs[@]}" --samplesLabel "${labels[@]}" \
        -R "$DNASE_DA/${set}_peaks.bed" \
        --beforeRegionStartLength 2000 --afterRegionStartLength 2000 \
        --skipZeros --numberOfProcessors "$SLURM_CPUS_PER_TASK" \
        -o "$CHIP_PLOTS/${set}_dars_chip_matrix.gz"

    plotHeatmap -m "$CHIP_PLOTS/${set}_dars_chip_matrix.gz" \
        -out "$CHIP_PLOTS/${set}_dars_chip_heatmap.png" \
        --colorMap Reds Reds Blues Blues Greens Greens \
        --regionsLabel "${set^} DNase peaks" \
        --refPointLabel centre --dpi 300

    plotProfile -m "$CHIP_PLOTS/${set}_dars_chip_matrix.gz" \
        -out "$CHIP_PLOTS/${set}_dars_chip_profile.png" --perGroup
done
