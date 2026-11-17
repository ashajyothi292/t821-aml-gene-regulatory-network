#!/bin/bash
# tests/run_tests.sh
#
# Runs the downstream steps on small SIMULATED data, so the code can be checked
# without the 2 GB sequencing files or the HPC cluster:
#   07_filter_peaks.sh -> 08_consensus_peaks.sh -> differential_accessibility.R
#   -> footprint_confidence_tiers.R -> 16_overlap_dars_with_chip.sh
#   -> deseq2_rnaseq.R -> gsea_preranked.py
# Alignment, peak calling, footprinting and HOMER need the real data and are not
# covered. The simulated numbers mean nothing biologically.
#
# Needs: bash, bedtools, R (DESeq2, GenomicRanges, dplyr, ggplot2, pheatmap,
#        org.Hs.eg.db), Python 3 (pandas, gseapy, matplotlib).
# Usage: bash tests/run_tests.sh   (from the repository root)

set -euo pipefail
REPO=$(pwd)
TMP=$(mktemp -d)
export PROJECT_DIR="$TMP"
module () { :; }          # HPC module commands do nothing here
export -f module
export SLURM_SUBMIT_DIR="$REPO"

echo "== simulating data in $TMP"
python3 "$REPO/tests/simulate_data.py" "$TMP"

echo "== peak filtering and consensus"
bash "$REPO/scripts/02_dnase_peaks/07_filter_peaks.sh"
bash "$REPO/scripts/02_dnase_peaks/08_consensus_peaks.sh"

echo "== stand-in for the HOMER annotation table"
python3 "$REPO/tests/simulate_data.py" "$TMP" --homer-table

echo "== differential accessibility"
Rscript "$REPO/scripts/03_differential_accessibility/differential_accessibility.R" "$TMP"

echo "== footprint tiers"
python3 "$REPO/tests/simulate_data.py" "$TMP" --footprints-and-chip
Rscript "$REPO/scripts/04_footprinting/footprint_confidence_tiers.R" "$TMP"

echo "== overlap with ChIP-seq"
bash "$REPO/scripts/05_chipseq_integration/16_overlap_dars_with_chip.sh"

echo "== RNA-seq"
Rscript "$REPO/scripts/06_rnaseq/deseq2_rnaseq.R" "$TMP"
python3 "$REPO/tests/simulate_data.py" "$TMP" --gmt
python3 "$REPO/scripts/06_rnaseq/gsea_preranked.py" "$TMP" "$TMP/reference/test_reactome.gmt"

echo "== checks"
fail=0
check () { if eval "$2"; then echo "  PASS  $1"; else echo "  FAIL  $1"; fail=1; fi; }
DA="$TMP/dnase/differential_accessibility"
check "consensus peaks written"            "[ -s $TMP/dnase/consensus/final_consensus_peaks.bed ]"
check "blacklisted summit removed"         "! grep -q '^chr1[[:space:]]5000[[:space:]]' $TMP/dnase/peaks/PBSC-1.height_filtered.bed"
check "non-standard chromosome removed"    "! grep -q chrUn $TMP/dnase/peaks/PBSC-1.height_filtered.bed"
check "gained regions found"               "[ \$(wc -l < $DA/gained_peaks.bed) -gt 50 ]"
check "lost regions found"                 "[ \$(wc -l < $DA/lost_peaks.bed) -gt 50 ]"
check "tier 1 gained regions written"      "[ -s $DA/tier1_gained.bed ]"
check "ChIP overlap table has 8 rows"      "[ \$(tail -n +2 $DA/dar_chip_overlap_counts.tsv | wc -l) -eq 8 ]"
check "regulatory elements written"        "[ -s $DA/gained_regulatory_elements.bed ]"
check "volcano plots drawn"                "[ -s $DA/volcano_differential_accessibility.png ] && [ -s $DA/volcano_footprint_confidence.png ]"
check "three GSEA rank files"              "[ \$(ls $TMP/rnaseq/deseq2/ranks_*.rnk | wc -l) -eq 3 ]"
check "RNA-seq plots drawn"                "[ -s $TMP/rnaseq/deseq2/pca.png ] && [ -s $TMP/rnaseq/deseq2/volcano.png ]"
check "FPKM table written"                 "[ -s $TMP/rnaseq/deseq2/fpkm.tsv ]"
check "GSEA results written"               "[ -s $TMP/rnaseq/gsea/reactome_significant.csv ]"

rm -rf "$TMP"
if [ $fail -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; exit 1; fi
