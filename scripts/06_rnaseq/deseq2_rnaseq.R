# deseq2_rnaseq.R
#
# Differential gene expression: t(8;21) AML patients vs healthy PBSC controls.
#
# Input:  featureCounts_gene_counts.txt (scripts/06_rnaseq/20_featurecounts.sh)
# Output: (in rnaseq/deseq2)
#         DESeq_RNA_all.csv              results for all six patients vs control
#         normalised_counts.csv
#         fpkm.tsv                       expression per gene, used to colour the networks
#         pca.png, volcano.png, sample_distances.png, top50_degs_heatmap.png
#         ranks_<comparison>.rnk          ranked gene lists for GSEA (gsea_preranked.py)
#
# Three comparisons are run, because the patients split into two clusters in
# the DEG / footprint heatmap:
#   all    all six t(8;21) patients vs control
#   core   patients without a co-occurring mutation or relapse (t821-1, -3, -4)
#   other  relapse or co-mutated samples (t821-1R, t821-ITDx2-5, t821-KIT-2)
#
# Usage:  Rscript scripts/06_rnaseq/deseq2_rnaseq.R <project_dir>

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) args[1] else "."
count_file <- file.path(project_dir, "rnaseq/counts/featureCounts_gene_counts.txt")
out_dir <- file.path(project_dir, "rnaseq/deseq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("PBSC-1", "PBSC-2", "t821-1", "t821-1R", "t821-3", "t821-4", "t821-ITDx2-5", "t821-KIT-2")
group <- c("Control", "Control", "core", "other", "core", "core", "other", "other")

# ---- counts ---------------------------------------------------------------------
fc <- read.delim(count_file, comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid
colnames(counts) <- samples            # featureCounts keeps the BAM order from the script

gene_info <- data.frame(
  ensembl = sub("\\..*$", "", fc$Geneid),
  chr     = sub(";.*$", "", fc$Chr),   # multi-exon genes list the chromosome once per exon
  row.names = fc$Geneid)
gene_info$symbol <- mapIds(org.Hs.eg.db, keys = gene_info$ensembl, keytype = "ENSEMBL",
                           column = "SYMBOL", multiVals = "first")

# ---- all patients vs control -----------------------------------------------------
coldata <- data.frame(row.names = samples,
                      condition = factor(ifelse(group == "Control", "Control", "t821"),
                                         levels = c("Control", "t821")),
                      group = group)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "t821", "Control"))

res_df <- data.frame(gene_id = rownames(res), symbol = gene_info[rownames(res), "symbol"],
                     as.data.frame(res), row.names = NULL)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(out_dir, "DESeq_RNA_all.csv"), row.names = FALSE)
write.csv(counts(dds, normalized = TRUE), file.path(out_dir, "normalised_counts.csv"))

# FPKM (fragments per kilobase of gene per million mapped fragments), the
# expression measure the gene regulatory network pipeline takes as input
gene_kb <- fc$Length / 1000
fpkm <- sweep(counts / gene_kb, 2, colSums(counts) / 1e6, "/")
fpkm_out <- data.frame(gene_id = rownames(counts), symbol = gene_info$symbol,
                       round(fpkm, 3), check.names = FALSE)
write.table(fpkm_out, file.path(out_dir, "fpkm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# DEGs: adjusted p < 0.1 and at least a two-fold change
res_df$status <- with(res_df, ifelse(!is.na(padj) & padj < 0.1 & log2FoldChange > 1, "Upregulated",
                              ifelse(!is.na(padj) & padj < 0.1 & log2FoldChange < -1, "Downregulated",
                                     "Not significant")))
cat("Genes tested:", sum(!is.na(res_df$padj)), "\n")
print(table(res_df$status))

# ---- plots ---------------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)

pca <- plotPCA(rld, intgroup = "condition") + theme_classic(base_size = 14) +
  ggtitle("PCA of AML and control samples")
ggsave(file.path(out_dir, "pca.png"), pca, width = 7, height = 5.5, dpi = 300)

volcano <- ggplot(res_df[!is.na(res_df$padj), ],
                  aes(log2FoldChange, -log10(padj), colour = status)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_colour_manual(values = c(Upregulated = "#E74C3C", Downregulated = "#2980B9",
                                 "Not significant" = "grey75")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(title = "Differential gene expression: AML vs control",
       subtitle = sprintf("Upregulated: %d  |  Downregulated: %d",
                          sum(res_df$status == "Upregulated"), sum(res_df$status == "Downregulated")),
       x = "Log2 fold change", y = "-Log10 adjusted p-value", colour = NULL)
ggsave(file.path(out_dir, "volcano.png"), volcano, width = 8, height = 6, dpi = 300)

dists <- as.matrix(dist(t(assay(rld))))
png(file.path(out_dir, "sample_distances.png"), width = 1800, height = 1600, res = 250)
pheatmap(dists, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Sample-to-sample distances (rlog)")
dev.off()

top50 <- head(res_df$gene_id[!is.na(res_df$padj)], 50)
mat <- assay(vsd)[top50, ]
rownames(mat) <- ifelse(is.na(gene_info[top50, "symbol"]), top50, gene_info[top50, "symbol"])
png(file.path(out_dir, "top50_degs_heatmap.png"), width = 1800, height = 2400, res = 250)
pheatmap(mat, scale = "row", annotation_col = coldata["condition"],
         main = "Top 50 genes by adjusted p-value (row z-score)")
dev.off()

# ---- ranked lists for GSEA ------------------------------------------------------------
# rank = -log10(adjusted p) x log2FC, after removing sex-chromosome, mitochondrial
# and ribosomal-protein genes (they dominate the ranking without saying anything
# about leukaemia biology)
write_ranks <- function(res, name) {
  df <- data.frame(gene_id = rownames(res), as.data.frame(res))
  df$symbol <- gene_info[df$gene_id, "symbol"]
  df$chr <- gene_info[df$gene_id, "chr"]
  df <- df[!is.na(df$padj) & !is.na(df$symbol) & df$padj > 0, ]
  df <- df[!df$chr %in% c("X", "Y", "MT") & !grepl("^(MT-|RPL|RPS)", df$symbol), ]
  df$rank <- -log10(df$padj) * df$log2FoldChange
  df <- df[order(abs(df$rank), decreasing = TRUE), ]
  df <- df[!duplicated(df$symbol), ]
  df <- df[order(df$rank, decreasing = TRUE), ]
  write.table(df[, c("symbol", "rank")], file.path(out_dir, sprintf("ranks_%s.rnk", name)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("ranks_%s.rnk: %d genes\n", name, nrow(df)))
}
write_ranks(res, "all")

for (sub in c("core", "other")) {
  keep <- group %in% c("Control", sub)
  cd <- data.frame(row.names = samples[keep],
                   condition = factor(ifelse(group[keep] == "Control", "Control", "t821"),
                                      levels = c("Control", "t821")))
  d <- DESeq(DESeqDataSetFromMatrix(counts[, keep], cd, design = ~ condition), quiet = TRUE)
  r <- results(d, contrast = c("condition", "t821", "Control"))
  write.csv(as.data.frame(r), file.path(out_dir, sprintf("DESeq_RNA_%s.csv", sub)))
  write_ranks(r, sub)
}
