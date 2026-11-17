#!/bin/bash
#SBATCH --job-name=footprint_motifs
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=bbdefault
#
# Which transcription factor motifs sit inside the regions that change?
#   1. motif enrichment in all gained / lost regions, against the consensus set
#   2. motif enrichment in the footprint-confirmed (Tier 1) gained / lost regions,
#      against the pooled t(8;21) or PBSC footprints
# Run after differential_accessibility.R and footprint_confidence_tiers.R
# (which writes the Tier 1 regions and the pooled footprint backgrounds).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load HOMER/4.11-foss-2020a

cd "$DNASE_DA"
mkdir -p homer_preparsed

run_homer () {   # <targets.bed> <background.bed> <output folder>
    findMotifsGenome.pl "$1" hg38 "$3" -size 200 -bg "$2" \
        -p 8 -preparsedDir homer_preparsed
}

run_homer gained_peaks.bed "$DNASE_CONSENSUS/final_consensus_peaks.bed" motifs_gained
run_homer lost_peaks.bed   "$DNASE_CONSENSUS/final_consensus_peaks.bed" motifs_lost
run_homer tier1_gained.bed t821_footprints.bed motifs_tier1_gained
run_homer tier1_lost.bed   PBSC_footprints.bed motifs_tier1_lost
