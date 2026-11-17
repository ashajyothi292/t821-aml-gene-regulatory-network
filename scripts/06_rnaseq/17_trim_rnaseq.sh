#!/bin/bash
#SBATCH --job-name=rna_trim
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=4
#SBATCH --qos=bbdefault
#SBATCH --array=0-7
#
# Trim paired-end RNA-seq reads with Trimmomatic: TruSeq adapters removed,
# 4-base sliding window with mean quality >= 15, reads shorter than 36 bp dropped.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2021b
module load Trimmomatic/0.39-Java-11

s=${RNA_SAMPLES[$SLURM_ARRAY_TASK_ID]}
mkdir -p "$RNA_TRIMMED"; cd "$RNA_TRIMMED"

java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" PE -threads 4 -phred33 \
    "$RNA_FASTQ/${s}_R1.fastq.gz" "$RNA_FASTQ/${s}_R2.fastq.gz" \
    "${s}_R1_trimmed_paired.fq.gz" "${s}_R1_trimmed_unpaired.fq.gz" \
    "${s}_R2_trimmed_paired.fq.gz" "${s}_R2_trimmed_unpaired.fq.gz" \
    ILLUMINACLIP:"$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
