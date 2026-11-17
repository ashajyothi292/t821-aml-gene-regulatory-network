#!/bin/bash
#SBATCH --job-name=dar_chip_overlap
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --qos=bbdefault
#
# Define high-confidence regulatory elements (Figure 11 of the report):
#   gained regions: does the region overlap RUNX1-ETO binding, H3K27ac (an
#                   active-enhancer mark), both, or neither?
#   lost regions:   the same with RUNX1 (displaced by RUNX1-ETO) and H3K27ac
# Regions bound by the factor AND marked by H3K27ac are kept as high-confidence
# elements; those that also contain a footprint (Tier 1) become the regulatory
# elements used to build the gene regulatory networks.
#
# Our counts: gained 11,739 -> 5,116 (44%) with RUNX1-ETO and H3K27ac;
#             lost   7,891 -> 3,374 (43%) with RUNX1 and H3K27ac.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2024a
module load BEDTools/2.31.1-GCC-13.3.0

cd "$DNASE_DA"
out=dar_chip_overlap_counts.tsv
printf "set\tcategory\tn_regions\n" > "$out"

classify () {   # <regions.bed> <factor peaks> <label> <factor name>
    local regions=$1 factor=$2 label=$3 fname=$4
    bedtools intersect -a "$regions" -b "$factor" -u > tmp_factor.bed
    bedtools intersect -a "$regions" -b "$CHIP_PEAKS/H3K27ac_peaks.bed" -u > tmp_k27.bed

    local total both factor_only k27_only neither
    total=$(wc -l < "$regions")
    both=$(bedtools intersect -a tmp_factor.bed -b "$CHIP_PEAKS/H3K27ac_peaks.bed" -u | tee "${label}_high_confidence.bed" | wc -l)
    factor_only=$(( $(wc -l < tmp_factor.bed) - both ))
    k27_only=$(( $(wc -l < tmp_k27.bed) - both ))
    neither=$(( total - both - factor_only - k27_only ))

    printf "%s\t%s only\t%d\n"            "$label" "$label" "$neither"  >> "$out"
    printf "%s\t%s + %s\t%d\n"            "$label" "$label" "$fname" "$factor_only" >> "$out"
    printf "%s\t%s + H3K27ac\t%d\n"       "$label" "$label" "$k27_only" >> "$out"
    printf "%s\tall three\t%d\n"          "$label" "$both"     >> "$out"
    rm tmp_factor.bed tmp_k27.bed
}

classify gained_peaks.bed "$CHIP_PEAKS/RUNX1-ETO_peaks.bed" gained RUNX1-ETO
classify lost_peaks.bed   "$CHIP_PEAKS/RUNX1_peaks.bed"     lost   RUNX1
cat "$out"

# final regulatory elements for the GRN: bound + H3K27ac + footprint confirmed
for label in gained lost; do
    bedtools intersect -a "${label}_high_confidence.bed" -b "tier1_${label}.bed" -u \
        > "${label}_regulatory_elements.bed"
    echo "${label} regulatory elements: $(wc -l < "${label}_regulatory_elements.bed")"
done
