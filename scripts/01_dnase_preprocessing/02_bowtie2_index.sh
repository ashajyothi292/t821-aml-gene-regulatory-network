#!/bin/bash
#SBATCH --job-name=bowtie2_index
#SBATCH --time=02:20:00
#SBATCH --mem=20G
#SBATCH --ntasks=10
#SBATCH --qos=bbdefault
#
# Build the Bowtie2 index for GRCh38 (needed once).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2023a
module load Bowtie2/2.5.4-GCC-12.3.0

mkdir -p "$(dirname "$BOWTIE2_INDEX")"
bowtie2-build --threads 10 "$REF_FASTA" "$BOWTIE2_INDEX"
