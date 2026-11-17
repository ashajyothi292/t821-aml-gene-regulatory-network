#!/bin/bash
#SBATCH --job-name=dedup
#SBATCH --time=03:20:00
#SBATCH --mem=30G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#SBATCH --array=0-6
#
# Remove PCR and optical duplicates with Picard MarkDuplicates, then index.
# Duplicates would otherwise inflate the accessibility signal at a site.
# OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 is the setting for patterned flow cells.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2023a
module load picard/3.4.0-Java-17

s=${DNASE_SAMPLES[$SLURM_ARRAY_TASK_ID]}
cd "$DNASE_BAM"

java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    -I "${s}.sorted.bam" \
    -O "${s}.dedup.bam" \
    -M "${s}.dup_metrics.txt" \
    -REMOVE_DUPLICATES true \
    -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

module purge; module load bluebear
module load bear-apps/2024a
module load SAMtools/1.21-GCC-13.3.0
samtools index "${s}.dedup.bam"
