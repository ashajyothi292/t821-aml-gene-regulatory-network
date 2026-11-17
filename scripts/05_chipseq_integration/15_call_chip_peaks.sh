#!/bin/bash
#SBATCH --job-name=chip_peaks
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=4
#SBATCH --qos=bbdefault
#
# Call binding sites from the ChIP-seq signal tracks, for the overlap analysis.
# GEO provides pileup bedGraphs rather than reads, so peaks are called from the
# tracks with MACS3:
#   bdgcmp -m ppois   Poisson -log10 p-value of each track against the Input
#   bdgpeakcall       regions above the cut-off (-c 1.0, i.e. p < 0.1),
#                     at least 500 bp long, gaps up to 1 kb joined
# Peaks from the 0Dox and 5Dox tracks of each factor are then merged.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2022b
module load MACS3/3.0.1-foss-2022b

mkdir -p "$CHIP_PEAKS"
for t in HA-RUNX1-ETO_0Dox HA-RUNX1-ETO_5Dox RUNX1_0Dox RUNX1_5Dox H3K27ac_0Dox H3K27ac_5Dox; do
    f=$(ls "$CHIP_DIR"/*"_${t}_treat_pileup.bed"[Gg]"raph" | head -n 1)
    macs3 bdgcmp -t "$f" -c "$CHIP_INPUT" -m ppois -o "$CHIP_PEAKS/${t}_ppois.bdg"
    macs3 bdgpeakcall -i "$CHIP_PEAKS/${t}_ppois.bdg" --cutoff 1.0 -l 500 -g 1000 \
        -o "$CHIP_PEAKS/${t}.narrowPeak"
done

module purge; module load bluebear
module load bear-apps/2024a
module load BEDTools/2.31.1-GCC-13.3.0

cd "$CHIP_PEAKS"
for factor in HA-RUNX1-ETO RUNX1 H3K27ac; do
    grep -v '^track' "${factor}_0Dox.narrowPeak" "${factor}_5Dox.narrowPeak" -h \
        | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > "${factor#HA-}_peaks.bed"
    echo "${factor#HA-}: $(wc -l < "${factor#HA-}_peaks.bed") merged peaks"
done
