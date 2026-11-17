#!/bin/bash
#SBATCH --job-name=filter_peaks
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=2
#SBATCH --qos=bbdefault
#
# Filter the MACS2 summits of every sample:
#   1. keep the standard chromosomes (chr1-22, X, Y)
#   2. drop summits inside ENCODE blacklisted regions
#   3. drop weak summits (MACS2 score column below 10)
# Writes <sample>.height_filtered.bed and a table of peak counts per sample.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2024a
module load BEDTools/2.31.1-GCC-13.3.0

cd "$DNASE_PEAKS"
: > peak_counts.tsv
for s in "${DNASE_SAMPLES[@]}"; do
    grep -E '^chr([0-9]+|X|Y)[[:space:]]' "${s}_summits.bed" \
        | bedtools intersect -a stdin -b "$BLACKLIST" -v \
        | awk 'BEGIN{OFS="\t"} $5 >= 10' > "${s}.height_filtered.bed"
    printf "%s\t%s\n" "$s" "$(wc -l < "${s}.height_filtered.bed")" >> peak_counts.tsv
done
cat peak_counts.tsv
