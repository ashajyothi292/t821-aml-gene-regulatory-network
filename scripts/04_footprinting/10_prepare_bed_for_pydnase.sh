#!/bin/bash
#SBATCH --job-name=prep_footprint_bed
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --qos=bbdefault
#
# Make the region files pyDNase expects: standard chromosomes only, three
# columns, Unix line endings, and each MACS2 summit widened by 50 bp per side
# so Wellington has room to look for a footprint.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

mkdir -p "$FOOTPRINT_DIR/regions"
for s in "${DNASE_SAMPLES[@]}"; do
    awk 'BEGIN{OFS="\t"}
         NF>=3 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
             start=$2-50; if (start<0) start=0
             print $1, start, $3+50
         }' "$DNASE_PEAKS/${s}_summits.bed" | tr -d '\r' \
        > "$FOOTPRINT_DIR/regions/${s}_summits.formatfixed.bed"
done
