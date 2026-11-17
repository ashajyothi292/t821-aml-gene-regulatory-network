#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --time=04:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=8
#SBATCH --qos=bbdefault
#
# Count read pairs per gene (Ensembl GRCh38 release 115) with featureCounts.
# -p: paired-end input   -s 0: unstranded library

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

mkdir -p "$RNA_COUNTS"
bams=(); for s in "${RNA_SAMPLES[@]}"; do bams+=("$RNA_BAM/${s}_Aligned.sortedByCoord.out.bam"); done

featureCounts -T 8 -p -s 0 \
    -a "$GTF" \
    -o "$RNA_COUNTS/featureCounts_gene_counts.txt" \
    "${bams[@]}"
