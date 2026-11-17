#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --time=02:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#SBATCH --array=0-7
#
# Align trimmed paired-end reads with STAR; output a coordinate-sorted BAM.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

s=${RNA_SAMPLES[$SLURM_ARRAY_TASK_ID]}
mkdir -p "$RNA_BAM"

STAR --genomeDir "$STAR_INDEX" --runThreadN 10 \
    --readFilesIn "$RNA_TRIMMED/${s}_R1_trimmed_paired.fq.gz" "$RNA_TRIMMED/${s}_R2_trimmed_paired.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$RNA_BAM/${s}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard
