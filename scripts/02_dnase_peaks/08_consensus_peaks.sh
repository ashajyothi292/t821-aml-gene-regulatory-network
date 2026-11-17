#!/bin/bash
#SBATCH --job-name=consensus_peaks
#SBATCH --time=01:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=2
#SBATCH --qos=bbdefault
#
# Build one consensus set of accessible regions across the cohort.
#   1. widen every filtered summit to a 401 bp window (summit +/- 200 bp)
#   2. merge all windows from all samples into a union of regions
#   3. keep a region if it is open in >= 3 of the 5 t(8;21) patients,
#      or in >= 1 of the 2 PBSC controls
#   4. combine the two lists -> final_consensus_peaks.bed
# In our data this gave 97,368 consensus regions.
#
# Also writes t821_union_peaks.bed (the patient DHS merged into one file),
# used as the reference regions for ChIP-seq signal (scripts 05_*).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

module purge; module load bluebear
module load bear-apps/2024a
module load BEDTools/2.31.1-GCC-13.3.0

mkdir -p "$DNASE_CONSENSUS"
cd "$DNASE_CONSENSUS"

for s in "${DNASE_SAMPLES[@]}"; do
    # summit position = start + summit offset; MACS2 summits files are 1 bp already
    awk 'BEGIN{OFS="\t"} {s=$2; print $1, (s>200 ? s-200 : 0), s+201}' \
        "$DNASE_PEAKS/${s}.height_filtered.bed" > "${s}.extended.bed"
done

cat ./*.extended.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > peak_union.bed
echo "union regions: $(wc -l < peak_union.bed)"

patient_beds=(); for s in "${DNASE_PATIENTS[@]}"; do patient_beds+=("${s}.extended.bed"); done
control_beds=(); for s in "${DNASE_CONTROLS[@]}"; do control_beds+=("${s}.extended.bed"); done

# -C reports, for each union region, how many files overlap it (one row per file)
bedtools intersect -a peak_union.bed -b "${patient_beds[@]}" -C \
    | awk 'BEGIN{OFS="\t"} $NF>0 {n[$1"\t"$2"\t"$3]++} END{for (k in n) if (n[k] >= 3) print k}' \
    | sort -k1,1 -k2,2n > t821_consensus.bed
bedtools intersect -a peak_union.bed -b "${control_beds[@]}" -C \
    | awk 'BEGIN{OFS="\t"} $NF>0 {n[$1"\t"$2"\t"$3]++} END{for (k in n) if (n[k] >= 1) print k}' \
    | sort -k1,1 -k2,2n > pbsc_consensus.bed

cat t821_consensus.bed pbsc_consensus.bed | sort -k1,1 -k2,2n | uniq > final_consensus_peaks.bed
echo "t(8;21) consensus: $(wc -l < t821_consensus.bed)"
echo "PBSC consensus:    $(wc -l < pbsc_consensus.bed)"
echo "final consensus:   $(wc -l < final_consensus_peaks.bed)"

# patient DHS merged into a single file, as suggested by our supervisor
cat "${patient_beds[@]}" | sort -k1,1 -k2,2n | bedtools merge -i stdin > t821_union_peaks.bed
