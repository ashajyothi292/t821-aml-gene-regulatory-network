#!/bin/bash
#SBATCH --job-name=dnase_trim
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=4
#SBATCH --qos=bbdefault
#SBATCH --array=0-6
#
# Adapter and quality trimming of single-end DNaseI-seq reads with Trimmomatic.
# One array task per sample. Submit from the repository root:
#   sbatch scripts/01_dnase_preprocessing/01_trim_dnase.sh
#
# Settings match the RNA-seq trimming (sliding window 4:15, minimum
# length 36), run here in single-end mode. FastQC was run before and after.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2021b
module load Trimmomatic/0.39-Java-11

s=${DNASE_SAMPLES[$SLURM_ARRAY_TASK_ID]}
mkdir -p "$DNASE_TRIMMED"

java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" SE -threads 4 -phred33 \
    "$DNASE_FASTQ/${s}.fastq.gz" \
    "$DNASE_TRIMMED/${s}_trimmed_SE.fastq.gz" \
    ILLUMINACLIP:"$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE.fa":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
