# Gene regulatory network construction

The networks were built with the gene regulatory network pipeline of
Coleman et al. (2023), *Cell Reports* 42(12):113568, which our supervisor's
group developed. Its scripts (`findMotifs.py`, `build_gene_regulatory_network.py`)
and the curated motif library belong to that group and are not included here.

## What goes in

| Input | Made by |
|---|---|
| Regulatory elements: `gained_regulatory_elements.bed`, `lost_regulatory_elements.bed` | `05_chipseq_integration/16_overlap_dars_with_chip.sh` |
| Footprint-confirmed regions for motif search: `tier1_gained.bed`, `tier1_lost.bed` | `04_footprinting/footprint_confidence_tiers.R` |
| Gene expression (FPKM per gene): `fpkm.tsv` | `06_rnaseq/deseq2_rnaseq.R` |
| Transcription factor motif library (position weight matrices) | Assi et al. (2019) |
| Transcription factor family annotation | pipeline resource |

## Steps

1. **Motif positions.** `findMotifs.py` scans the footprint-confirmed regions
   (hg38) with the motif library and writes one BED file of motif positions per
   transcription factor.
2. **Network.** `build_gene_regulatory_network.py` links each motif in a
   regulatory element to the gene that element is assigned to, keeps genes that
   are expressed, and writes transcription factor to target edges, with an
   interaction count per edge. Run once for gained and once for lost elements.
3. **Visualisation.** The `.cyjs` output was opened in Cytoscape. Node size =
   number of downstream factors (out-degree), node colour = expression (FPKM),
   edge colour = interaction count; factors grouped by family.

## Result

| Network | Nodes | Edges |
|---|---|---|
| Gained elements | 63 | 214 |
| Lost elements | 61 | 135 |

Factors such as MYC, FOS and ELF1 appear in both networks with different
connections, which points to re-wiring of existing factors at newly opened
regions rather than a completely different set of regulators.
