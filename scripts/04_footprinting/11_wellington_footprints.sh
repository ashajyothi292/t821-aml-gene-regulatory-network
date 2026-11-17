#!/bin/bash
#SBATCH --job-name=wellington
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=bbdefault
#SBATCH --array=0-6
#
# DNase I footprinting with the Wellington algorithm (pyDNase).
# A footprint is a short stretch inside an open region that DNase I could not
# cut because a transcription factor was sitting on the DNA.
# -fdrlimit -4 keeps footprints with a p-value below 1e-4; the
# WellingtonFootprints.FDR.0.01.bed output is used downstream.
# Also writes forward / reverse cut-count tracks for plotting.
#
# Footprints found per sample in our data: PBSC-1 3,426; PBSC-2 670;
# t821-1 4,055; t821-1rel 15,427; t821-3 8,949; t821-4 5,123; t821-KIT 4,478.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load pyDNase/0.3.0-foss-2019b-Python-3.7.4
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

s=${DNASE_SAMPLES[$SLURM_ARRAY_TASK_ID]}
regions="$FOOTPRINT_DIR/regions/${s}_summits.formatfixed.bed"
bam="$DNASE_BAM/${s}.dedup.bam"
out="$FOOTPRINT_DIR/${s}"
mkdir -p "$out"; cd "$out"

wellington_footprints.py "$regions" "$bam" . -fdrlimit -4
dnase_wig_tracks.py "$regions" "$bam" "${s}_cuts_fw.wig" "${s}_cuts_rev.wig"

# average cut profile around the footprints (a QC plot: a dip in the middle
# means the footprints are real)
dnase_average_profile.py -w 200 \
    ./*WellingtonFootprints.FDR.0.01.bed "$bam" "${s}_average_profile"
