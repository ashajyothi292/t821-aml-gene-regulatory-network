# footprint_confidence_tiers.R
#
# Not every region that changes accessibility is actually bound by a
# transcription factor. Here each differentially accessible region (DAR) is
# checked for a Wellington footprint:
#   Tier 1  DAR that contains a footprint (evidence of factor binding)
#   Tier 2  DAR without a footprint
#   Tier 3  unchanged region
# Gained regions are checked against t(8;21) footprints, lost regions against
# PBSC footprints. Tier 1 regions are the input for motif analysis and the
# gene regulatory network.
#
# In our data 4,835 of 11,739 gained (41.2%) and 3,274 of 7,891 lost (41.5%)
# regions were Tier 1.
#
# Usage: Rscript scripts/04_footprinting/footprint_confidence_tiers.R <project_dir>
# Needs: differential_accessibility.R and 11_wellington_footprints.sh to have run.
# Also writes the pooled footprints of each condition (t821_footprints.bed,
# PBSC_footprints.bed), used as motif backgrounds in 12_footprint_motifs.sh.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) args[1] else "."
da_dir <- file.path(project_dir, "dnase/differential_accessibility")

res <- read.csv(file.path(da_dir, "DESeq2_results_annotated.csv"))

# ---- pool the Wellington footprints of each condition ------------------------------
fp_dir <- file.path(project_dir, "dnase/footprinting")
pool_footprints <- function(prefix) {
  files <- Sys.glob(file.path(fp_dir, paste0(prefix, "*"), "*WellingtonFootprints.FDR.0.01.bed"))
  if (!length(files)) stop("No Wellington footprint files found for ", prefix, " in ", fp_dir)
  fp <- do.call(rbind, lapply(files, function(f) read.table(f)[, 1:3]))
  colnames(fp) <- c("Chr", "Start", "End")
  fp <- fp[order(fp$Chr, fp$Start), ]
  write.table(fp, file.path(da_dir, paste0(sub("-$", "", prefix), "_footprints.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  makeGRangesFromDataFrame(fp, seqnames.field = "Chr", start.field = "Start", end.field = "End")
}
t821_fp <- pool_footprints("t821-")
pbsc_fp <- pool_footprints("PBSC-")
cat("Pooled footprints: t(8;21)", length(t821_fp), "| PBSC", length(pbsc_fp), "\n")

res_gr <- makeGRangesFromDataFrame(res, seqnames.field = "Chr", start.field = "Start",
                                   end.field = "End", keep.extra.columns = TRUE)

res$has_t821_footprint <- overlapsAny(res_gr, t821_fp)
res$has_pbsc_footprint <- overlapsAny(res_gr, pbsc_fp)

res <- res %>%
  mutate(Confidence = case_when(
    Category == "Gained (t821)" &  has_t821_footprint ~ "Tier1_Gained_Footprint",
    Category == "Lost (t821)"   &  has_pbsc_footprint ~ "Tier1_Lost_Footprint",
    Category == "Gained (t821)"                       ~ "Tier2_Gained_Only",
    Category == "Lost (t821)"                         ~ "Tier2_Lost_Only",
    TRUE                                              ~ "Tier3_Unchanged"))

tier_summary <- count(res, Confidence)
print(tier_summary)
cat(sprintf("Gained regions with a footprint: %.1f%%\n",
            100 * mean(res$has_t821_footprint[res$Category == "Gained (t821)"])))
cat(sprintf("Lost regions with a footprint:   %.1f%%\n",
            100 * mean(res$has_pbsc_footprint[res$Category == "Lost (t821)"])))

write.csv(res, file.path(da_dir, "DESeq2_results_with_footprints.csv"), row.names = FALSE)

# ---- Tier 1 regions as BED, for HOMER and the GRN pipeline ------------------------
for (tier in c("gained", "lost")) {
  label <- if (tier == "gained") "Tier1_Gained_Footprint" else "Tier1_Lost_Footprint"
  bed <- res %>%
    filter(Confidence == label) %>%
    arrange(Chr, Start) %>%
    transmute(Chr, Start, End, Name = paste(Chr, Start, End, sep = "_"), Score = 0, Strand = "+")
  write.table(bed, file.path(da_dir, sprintf("tier1_%s.bed", tier)), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ---- genes with the most Tier 1 regions -------------------------------------------
top_genes <- function(label) {
  res %>%
    filter(Confidence == label, !is.na(GeneName), GeneName != "") %>%
    group_by(GeneName) %>%
    summarise(n_regions = n(), mean_log2FC = mean(log2FoldChange), best_padj = min(padj),
              .groups = "drop") %>%
    arrange(desc(n_regions), best_padj)
}
write.csv(top_genes("Tier1_Gained_Footprint"), file.path(da_dir, "tier1_gained_genes.csv"), row.names = FALSE)
write.csv(top_genes("Tier1_Lost_Footprint"),   file.path(da_dir, "tier1_lost_genes.csv"),   row.names = FALSE)

# ---- volcano plot coloured by tier ------------------------------------------------
volcano <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = Confidence)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = c(Tier1_Gained_Footprint = "#C0392B", Tier1_Lost_Footprint = "#1A5276",
                                 Tier2_Gained_Only = "#E74C3C", Tier2_Lost_Only = "#2980B9",
                                 Tier3_Unchanged = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 14) +
  labs(title = "Differential accessibility with footprint confidence",
       subtitle = "Dark = footprint confirmed  |  Light = accessibility change only",
       x = "Log2 fold change", y = "-Log10 adjusted p-value", colour = "Confidence tier")
ggsave(file.path(da_dir, "volcano_footprint_confidence.png"), volcano, width = 10, height = 7, dpi = 300)
