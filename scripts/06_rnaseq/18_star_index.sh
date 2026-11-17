#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=03:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=20
#SBATCH --qos=bbdefault
#
# Build the STAR genome index for GRCh38 with Ensembl gene models (needed once).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

mkdir -p "$STAR_INDEX"
gtf_plain="$STAR_INDEX/annotation.gtf"
zcat "$GTF" > "$gtf_plain"

STAR --runMode genomeGenerate --runThreadN 20 \
    --genomeDir "$STAR_INDEX" \
    --genomeFastaFiles "$REF_FASTA" \
    --sjdbGTFfile "$gtf_plain"
