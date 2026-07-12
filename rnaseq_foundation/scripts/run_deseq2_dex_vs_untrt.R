suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

message("Task 10: DESeq2 differential expression analysis")

txi_path <- "results/tximport/txi_salmon_gene_level.rds"
samples_path <- "metadata/samples.csv"
annotation_path <- "results/tximport/gene_annotation.tsv"

if (!file.exists(txi_path)) stop("Missing tximport RDS: ", txi_path)
if (!file.exists(samples_path)) stop("Missing samples.csv: ", samples_path)
if (!file.exists(annotation_path)) stop("Missing gene annotation: ", annotation_path)

txi <- readRDS(txi_path)
samples <- read.csv(samples_path, stringsAsFactors = FALSE)

required_cols <- c("sample_id", "cell_line", "dex", "condition")
missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("metadata/samples.csv is missing columns: ", paste(missing_cols, collapse = ", "))
}

# Match sample metadata to tximport count matrix columns
samples <- samples[match(colnames(txi$counts), samples$sample_id), ]

if (!all(samples$sample_id == colnames(txi$counts))) {
  stop("Sample order mismatch between tximport counts and metadata/samples.csv")
}

rownames(samples) <- samples$sample_id

# Important paired design:
# dex effect is tested while controlling for cell_line.
samples$cell_line <- factor(samples$cell_line)
samples$dex <- factor(samples$dex, levels = c("untrt", "trt"))

message("Sample metadata used:")
print(samples)

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

write.table(
  samples,
  file = "results/deseq2/sample_metadata_used.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Creating DESeqDataSet from tximport...")
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ cell_line + dex
)

message("Genes before filtering: ", nrow(dds))

# Basic low-count filtering:
# keep genes with count >= 10 in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

message("Genes after filtering: ", nrow(dds))

message("Running DESeq2...")
dds <- DESeq(dds)

message("Available result names:")
print(resultsNames(dds))

res <- results(dds, contrast = c("dex", "trt", "untrt"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- res_df[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

annotation <- read.delim(annotation_path, stringsAsFactors = FALSE)
annotation <- unique(annotation[, c("gene_id", "gene_name", "gene_type")])

# A gene_id may appear with multiple gene_type values because the annotation
# was derived from transcript-level GENCODE records.
# Collapse annotation to one row per gene before merging with DESeq2 results.
collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  paste(x, collapse = ";")
}

annotation_gene <- aggregate(
  cbind(gene_name, gene_type) ~ gene_id,
  data = annotation,
  FUN = collapse_unique
)

res_df <- merge(annotation_gene, res_df, by = "gene_id", all.y = TRUE, sort = FALSE)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.table(
  res_df,
  file = "results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_padj <- subset(res_df, !is.na(padj) & padj < 0.05)
sig_padj_lfc1 <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

write.table(
  sig_padj,
  file = "results/deseq2/significant_genes_padj_0.05.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sig_padj_lfc1,
  file = "results/deseq2/significant_genes_padj_0.05_abs_log2fc_1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- data.frame(
  gene_id = rownames(norm_counts),
  norm_counts,
  check.names = FALSE
)

write.table(
  norm_counts_df,
  file = "results/deseq2/normalized_counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(dds, file = "results/deseq2/dds_deseq2.rds")
saveRDS(res, file = "results/deseq2/res_dex_trt_vs_untrt.rds")

summary_df <- data.frame(
  item = c(
    "samples",
    "genes_before_filtering",
    "genes_after_filtering",
    "significant_padj_0.05",
    "significant_padj_0.05_up",
    "significant_padj_0.05_down",
    "significant_padj_0.05_abs_log2fc_1"
  ),
  value = c(
    nrow(samples),
    nrow(txi$counts),
    nrow(dds),
    nrow(sig_padj),
    sum(sig_padj$log2FoldChange > 0, na.rm = TRUE),
    sum(sig_padj$log2FoldChange < 0, na.rm = TRUE),
    nrow(sig_padj_lfc1)
  )
)

write.table(
  summary_df,
  file = "results/deseq2/task10_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Creating PCA and QC plots...")

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("dex", "cell_line"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

write.table(
  pca_data,
  file = "results/deseq2/pca_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf("results/figures/task10_pca_plot.pdf", width = 7, height = 5)
print(
  ggplot(pca_data, aes(x = PC1, y = PC2, color = dex, shape = cell_line, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle("PCA of airway RNA-seq samples") +
    theme_bw()
)
dev.off()

pdf("results/figures/task10_ma_plot.pdf", width = 7, height = 5)
plotMA(res, ylim = c(-6, 6), main = "DESeq2 MA plot: dex treated vs untreated")
dev.off()

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- samples$sample_id
colnames(sample_dist_matrix) <- samples$sample_id

write.table(
  sample_dist_matrix,
  file = "results/deseq2/sample_distance_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

annotation_col <- data.frame(
  dex = samples$dex,
  cell_line = samples$cell_line
)
rownames(annotation_col) <- samples$sample_id

pdf("results/figures/task10_sample_distance_heatmap.pdf", width = 7, height = 6)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annotation_col,
  main = "Sample distance heatmap"
)
dev.off()

message("Task 10 finished.")
message("Main outputs:")
message("  results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv")
message("  results/deseq2/significant_genes_padj_0.05.tsv")
message("  results/deseq2/task10_summary.tsv")
message("  results/figures/task10_pca_plot.pdf")
message("  results/figures/task10_ma_plot.pdf")
message("  results/figures/task10_sample_distance_heatmap.pdf")
