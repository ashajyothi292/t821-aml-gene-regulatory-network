# differential_accessibility.R
#
# Which DNase I hypersensitive sites (DHS) gain or lose accessibility in
# t(8;21) AML patients compared with healthy PBSC controls?
#
# Input:  peak_signal_annotated.txt from HOMER (scripts/02_dnase_peaks/09_*):
#         one row per consensus region, HOMER annotation columns, then one
#         pileup-signal column per sample
# Output: (in the differential_accessibility folder)
#         DESeq2_results_annotated.csv   every region with log2FC, padj and annotation
#         gained_peaks.csv / lost_peaks.csv and .bed
#         volcano_differential_accessibility.png
#         genomic_annotation_gained_lost.png
#
# Usage:  Rscript scripts/03_differential_accessibility/differential_accessibility.R <project_dir>

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) args[1] else "."
in_file <- file.path(project_dir, "dnase/consensus/peak_signal_annotated.txt")
out_dir <- file.path(project_dir, "dnase/differential_accessibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("PBSC-1", "PBSC-2", "t821-1", "t821-1rel", "t821-3", "t821-4", "t821-KIT")
condition <- factor(ifelse(grepl("^PBSC", samples), "PBSC", "t821"), levels = c("PBSC", "t821"))

# ---- read the HOMER table ------------------------------------------------------
homer <- read.delim(in_file, check.names = FALSE, quote = "")
colnames(homer)[1] <- "PeakID"

# the first 19 columns are HOMER's annotation; the rest are the bedGraph signals,
# in the order the files were given to annotatePeaks.pl
annot <- homer[, 1:19]
colnames(annot) <- c("PeakID", "Chr", "Start", "End", "Strand", "PeakScore", "FocusRatio",
                     "Annotation", "DetailedAnnotation", "DistToTSS", "NearestPromoterID",
                     "EntrezID", "NearestUnigene", "NearestRefseq", "NearestEnsembl",
                     "GeneName", "GeneAlias", "GeneDescription", "GeneType")
signal <- as.matrix(homer[, 20:(19 + length(samples))])
colnames(signal) <- samples
rownames(signal) <- annot$PeakID
signal[is.na(signal)] <- 0

cat("Consensus regions:", nrow(signal), "\n")

# ---- DESeq2 ---------------------------------------------------------------------
# DESeq2 works on integer counts, so the pileup signal is rounded. Its
# median-of-ratios size factors then normalise for sequencing depth.
counts <- round(signal)
storage.mode(counts) <- "integer"

coldata <- data.frame(row.names = samples, condition = condition)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "t821", "PBSC"))   # positive = more open in AML

res_df <- data.frame(PeakID = rownames(res), as.data.frame(res), row.names = NULL)

res_annotated <- res_df %>%
  inner_join(annot, by = "PeakID") %>%
  mutate(
    DistToTSS = as.numeric(DistToTSS),
    RegionType = if_else(abs(DistToTSS) < 1500, "Promoter-proximal", "Distal"),
    # HOMER writes e.g. "intron (NM_001, intron 2 of 5)"; keep the first word
    Feature = sub(" .*", "", Annotation),
    Feature = recode(Feature, "promoter-TSS" = "Promoter", "intron" = "Intron",
                     "Intergenic" = "Intergenic", "exon" = "Exon", "TTS" = "TTS",
                     .default = "Other"),
    Category = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Gained (t821)",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Lost (t821)",
      TRUE ~ "Not significant"))

write.csv(res_annotated, file.path(out_dir, "DESeq2_results_annotated.csv"), row.names = FALSE)

gained <- filter(res_annotated, Category == "Gained (t821)")
lost   <- filter(res_annotated, Category == "Lost (t821)")
write.csv(gained, file.path(out_dir, "gained_peaks.csv"), row.names = FALSE)
write.csv(lost,   file.path(out_dir, "lost_peaks.csv"),   row.names = FALSE)

write_bed <- function(df, path) {
  write.table(df[order(df$Chr, df$Start), c("Chr", "Start", "End", "PeakID")], path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
write_bed(gained, file.path(out_dir, "gained_peaks.bed"))
write_bed(lost,   file.path(out_dir, "lost_peaks.bed"))

cat("Gained in t(8;21):", nrow(gained), "\n")
cat("Lost in t(8;21):  ", nrow(lost), "\n")
print(table(res_annotated$RegionType))

# ---- volcano plot ------------------------------------------------------------------
volcano <- ggplot(res_annotated, aes(log2FoldChange, -log10(padj), colour = Category)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = c("Gained (t821)" = "#E74C3C", "Lost (t821)" = "#2980B9",
                                 "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 14) +
  labs(title = "Differential chromatin accessibility: t(8;21) vs PBSC",
       subtitle = sprintf("Gained: %d  |  Lost: %d", nrow(gained), nrow(lost)),
       x = "Log2 fold change", y = "-Log10 adjusted p-value", colour = NULL)
ggsave(file.path(out_dir, "volcano_differential_accessibility.png"), volcano,
       width = 9, height = 6.5, dpi = 300)

# ---- where do the changed regions sit? ----------------------------------------------
feature_pct <- bind_rows(gained, lost) %>%
  filter(Feature != "Other") %>%
  count(Category, Feature) %>%
  group_by(Category) %>%
  mutate(Percent = 100 * n / sum(n))

annotation_plot <- ggplot(feature_pct, aes(Feature, Percent, fill = Category)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  scale_fill_manual(values = c("Gained (t821)" = "#E74C3C", "Lost (t821)" = "#2980B9")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genomic annotation: gained vs lost regions",
       x = "Genomic feature", y = "Percentage of regions (%)", fill = NULL)
ggsave(file.path(out_dir, "genomic_annotation_gained_lost.png"), annotation_plot,
       width = 8, height = 6, dpi = 300)
