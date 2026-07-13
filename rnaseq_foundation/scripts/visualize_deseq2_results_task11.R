suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

message("Task 11: visualize and interpret DESeq2 results")

res_path <- "results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv"
dds_path <- "results/deseq2/dds_deseq2.rds"
samples_path <- "metadata/samples.csv"
norm_counts_path <- "results/deseq2/normalized_counts.tsv"

if (!file.exists(res_path)) stop("Missing DESeq2 result table: ", res_path)
if (!file.exists(dds_path)) stop("Missing DESeq2 dds object: ", dds_path)
if (!file.exists(samples_path)) stop("Missing samples.csv: ", samples_path)
if (!file.exists(norm_counts_path)) stop("Missing normalized counts: ", norm_counts_path)

dir.create("results/task11", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

res <- read.delim(res_path, stringsAsFactors = FALSE)
dds <- readRDS(dds_path)
samples <- read.csv(samples_path, stringsAsFactors = FALSE)
norm_counts <- read.delim(norm_counts_path, stringsAsFactors = FALSE, check.names = FALSE)

samples$dex <- factor(samples$dex, levels = c("untrt", "trt"))
samples$cell_line <- factor(samples$cell_line)
rownames(samples) <- samples$sample_id

message("Loaded DESeq2 results: ", nrow(res), " genes")

# -----------------------------
# 1. Volcano plot
# -----------------------------

res$padj_plot <- res$padj
res$padj_plot[is.na(res$padj_plot)] <- 1
res$padj_plot[res$padj_plot == 0] <- .Machine$double.xmin

res$neg_log10_padj <- -log10(res$padj_plot)

res$direction <- "Not significant"
res$direction[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 1] <- "Up"
res$direction[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange <= -1] <- "Down"

volcano_summary <- as.data.frame(table(res$direction))
colnames(volcano_summary) <- c("category", "count")

write.table(
  volcano_summary,
  file = "results/task11/volcano_category_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

label_genes <- res[!is.na(res$padj), ]
label_genes <- label_genes[order(label_genes$padj), ]
label_genes <- head(label_genes, 12)

pdf("results/figures/task11_volcano_plot.pdf", width = 7, height = 6)
print(
  ggplot(res, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text(
      data = label_genes,
      aes(label = gene_name),
      size = 3,
      vjust = -0.5,
      check_overlap = TRUE
    ) +
    xlab("log2 fold change: dex treated vs untreated") +
    ylab("-log10 adjusted p-value") +
    ggtitle("Volcano plot: dexamethasone treated vs untreated") +
    theme_bw()
)
dev.off()

# -----------------------------
# 2. Top genes table
# -----------------------------

candidate <- subset(
  res,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1 &
    baseMean >= 100
)

candidate$direction <- ifelse(candidate$log2FoldChange > 0, "up", "down")
candidate$linear_fold_change <- ifelse(
  candidate$log2FoldChange >= 0,
  2 ^ candidate$log2FoldChange,
  -1 * (2 ^ abs(candidate$log2FoldChange))
)

candidate <- candidate[order(candidate$padj), ]

write.table(
  candidate,
  file = "results/task11/candidate_genes_padj_0.05_abs_log2fc_1_basemean_100.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  head(candidate, 50),
  file = "results/task11/top50_candidate_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 3. Top 50 heatmap
# -----------------------------

top50 <- head(candidate$gene_id, 50)

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

top50 <- top50[top50 %in% rownames(vsd_mat)]
heat_mat <- vsd_mat[top50, , drop = FALSE]

# Z-score by gene across samples
heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z[is.na(heat_mat_z)] <- 0

gene_labels <- res$gene_name[match(rownames(heat_mat_z), res$gene_id)]
gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(heat_mat_z)[is.na(gene_labels) | gene_labels == ""]
gene_labels <- paste0(gene_labels, " | ", rownames(heat_mat_z))

annotation_col <- samples[colnames(heat_mat_z), c("dex", "cell_line"), drop = FALSE]

pdf("results/figures/task11_top50_candidate_gene_heatmap.pdf", width = 8, height = 10)
pheatmap(
  heat_mat_z,
  labels_row = gene_labels,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  fontsize_row = 6,
  main = "Top 50 DE genes after dexamethasone treatment"
)
dev.off()

# -----------------------------
# 4. Selected gene expression plots
# -----------------------------

selected_gene_names <- c(
  "ZBTB16", "DUSP1", "PER1", "VCAM1",
  "SAMHD1", "GPX3", "SERPINA3", "ADAMTS1"
)

selected_res <- res[res$gene_name %in% selected_gene_names, ]
selected_res <- selected_res[order(match(selected_res$gene_name, selected_gene_names)), ]

write.table(
  selected_res,
  file = "results/task11/selected_gene_deseq2_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

selected_gene_ids <- selected_res$gene_id
selected_gene_ids <- selected_gene_ids[selected_gene_ids %in% norm_counts$gene_id]

plot_df_list <- list()

for (gid in selected_gene_ids) {
  row <- norm_counts[norm_counts$gene_id == gid, ]
  gene_name <- res$gene_name[match(gid, res$gene_id)]

  values <- as.numeric(row[1, samples$sample_id])
  df <- data.frame(
    gene_id = gid,
    gene_name = gene_name,
    sample_id = samples$sample_id,
    cell_line = samples$cell_line,
    dex = samples$dex,
    normalized_count = values
  )
  plot_df_list[[gid]] <- df
}

plot_df <- do.call(rbind, plot_df_list)

write.table(
  plot_df,
  file = "results/task11/selected_gene_normalized_counts_long.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf("results/figures/task11_selected_gene_expression.pdf", width = 10, height = 7)
print(
  ggplot(plot_df, aes(x = dex, y = log10(normalized_count + 1), group = cell_line)) +
    geom_line(alpha = 0.7) +
    geom_point(aes(shape = cell_line), size = 2.5) +
    facet_wrap(~ gene_name, scales = "free_y") +
    xlab("Treatment") +
    ylab("log10 normalized count + 1") +
    ggtitle("Selected dexamethasone-responsive genes") +
    theme_bw()
)
dev.off()

# -----------------------------
# 5. Summary
# -----------------------------

summary_df <- data.frame(
  item = c(
    "genes_in_result_table",
    "padj_0.05",
    "padj_0.05_up",
    "padj_0.05_down",
    "padj_0.05_abs_log2fc_1",
    "candidate_padj_0.05_abs_log2fc_1_basemean_100",
    "selected_genes_plotted"
  ),
  value = c(
    nrow(res),
    sum(!is.na(res$padj) & res$padj < 0.05),
    sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 0),
    sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < 0),
    sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1),
    nrow(candidate),
    length(unique(plot_df$gene_name))
  )
)

write.table(
  summary_df,
  file = "results/task11/task11_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Task 11 finished.")
message("Main outputs:")
message("  results/figures/task11_volcano_plot.pdf")
message("  results/figures/task11_top50_candidate_gene_heatmap.pdf")
message("  results/figures/task11_selected_gene_expression.pdf")
message("  results/task11/candidate_genes_padj_0.05_abs_log2fc_1_basemean_100.tsv")
message("  results/task11/task11_summary.tsv")
