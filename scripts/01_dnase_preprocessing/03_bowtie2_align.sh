#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --time=07:20:00
#SBATCH --mem=100G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#SBATCH --array=0-6
#
# Align trimmed single-end DNaseI-seq reads to GRCh38 with Bowtie2.
# --very-sensitive-local soft-clips read ends, which suits short DNase fragments;
# --no-unal keeps unaligned reads out of the SAM file.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2023a
module load Bowtie2/2.5.4-GCC-12.3.0

s=${DNASE_SAMPLES[$SLURM_ARRAY_TASK_ID]}
mkdir -p "$DNASE_BAM"

bowtie2 --no-unal --very-sensitive-local -p 6 \
    -x "$BOWTIE2_INDEX" \
    -U "$DNASE_TRIMMED/${s}_trimmed_SE.fastq.gz" \
    -S "$DNASE_BAM/${s}.sam"
