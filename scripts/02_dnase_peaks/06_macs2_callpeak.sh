#!/bin/bash
#SBATCH --job-name=macs2_dnase
#SBATCH --time=03:20:00
#SBATCH --mem=50G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#SBATCH --array=0-6
#
# Call DNase I hypersensitive sites (DHS) for each sample with MACS2.
# --nomodel:     DNase cuts are not ChIP fragments, so no fragment-size model
# --keep-dup all: duplicates were already removed by Picard
# -q 0.0005:     strict cut-off for high-confidence sites
# --call-summits / -B: report summits and write pileup bedGraphs, which are
#                     reused for quantification (script 09)

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2022a
module load MACS2/2.2.9.1-foss-2022a

s=${DNASE_SAMPLES[$SLURM_ARRAY_TASK_ID]}
mkdir -p "$DNASE_PEAKS"

macs2 callpeak -t "$DNASE_BAM/${s}.dedup.bam" -f BAM -g hs \
    --keep-dup all --nomodel -q 0.0005 --call-summits \
    -B --trackline -n "$s" --outdir "$DNASE_PEAKS"
