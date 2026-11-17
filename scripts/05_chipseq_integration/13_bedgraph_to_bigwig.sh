#!/bin/bash
#SBATCH --job-name=chip_bigwig
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#
# Convert the published ChIP-seq signal tracks (GEO GSE137670) to bigWig.
#
# The data come from Nafria et al. (2020): a doxycycline-inducible RUNX1-ETO
# model of early human myeloid precursors, before (0Dox) and after 5 days of
# induction (5Dox), for RUNX1-ETO, RUNX1 and H3K27ac (plus H3K4me3, LMO2, LDB1
# and an Input control). Download the *_treat_pileup.bedgraph files from GEO
# into $CHIP_DIR first.
#
# Each file is cleaned (UCSC track/browser lines removed), sorted and converted
# with bedGraphToBigWig. bigWig is what deepTools reads in the next step.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2023a
module load Kent_tools/478-GCC-12.3.0

mkdir -p "$CHIP_BW"
shopt -s nullglob
for f in "$CHIP_DIR"/*_treat_pileup.bedgraph "$CHIP_DIR"/*_treat_pileup.bedGraph; do
    name=$(basename "$f" | sed -E 's/_treat_pileup\.bed[Gg]raph$//')
    echo "converting $name"
    # bedGraphToBigWig needs plain byte-order (C locale) sorting
    grep -v -E '^(track|browser|type)' "$f" \
        | LC_ALL=C sort -k1,1 -k2,2n > "$CHIP_BW/${name}.sorted.bedgraph"
    bedGraphToBigWig "$CHIP_BW/${name}.sorted.bedgraph" "$CHROM_SIZES" "$CHIP_BW/${name}.bw"
    rm "$CHIP_BW/${name}.sorted.bedgraph"
done
ls -lh "$CHIP_BW"
