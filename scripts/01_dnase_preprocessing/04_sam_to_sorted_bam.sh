#!/bin/bash
#SBATCH --job-name=sam_to_bam
#SBATCH --time=03:20:00
#SBATCH --mem=30G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#SBATCH --array=0-6
#
# Convert SAM to BAM, sort by coordinate and report mapping statistics.
# The large SAM file is deleted once the sorted BAM exists.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2024a
module load SAMtools/1.21-GCC-13.3.0

s=${DNASE_SAMPLES[$SLURM_ARRAY_TASK_ID]}
cd "$DNASE_BAM"

samtools view -@ 6 -b -o "${s}.bam" "${s}.sam"
samtools sort -@ 6 -o "${s}.sorted.bam" "${s}.bam"
samtools flagstat "${s}.sorted.bam" > "${s}.flagstat.txt"

rm "${s}.sam" "${s}.bam"
echo "${s} done"
