#!/bin/bash
#SBATCH --job-name=homer_annotate
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=4
#SBATCH --qos=bbdefault
#
# Annotate the consensus regions with HOMER (nearest gene, distance to TSS,
# promoter / intron / intergenic ...) and, in the same pass, read the MACS2
# pileup signal of every sample over each region (one -bedGraph per sample).
# The output table is the input for differential accessibility in R.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load HOMER/4.11-foss-2020a

cd "$DNASE_CONSENSUS"

# HOMER needs a unique ID per region
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak_"NR,".","+"}' final_consensus_peaks.bed \
    > final_consensus_peaks.homer.bed

bedgraphs=()
for s in "${DNASE_SAMPLES[@]}"; do bedgraphs+=("$DNASE_PEAKS/${s}_treat_pileup.bdg"); done

annotatePeaks.pl final_consensus_peaks.homer.bed hg38 -size given \
    -bedGraph "${bedgraphs[@]}" > peak_signal_annotated.txt

echo "annotated regions: $(($(wc -l < peak_signal_annotated.txt) - 1))"
