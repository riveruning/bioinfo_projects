suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

message("Task 12: GO enrichment and GSEA")

res_path <- "results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv"

if (!file.exists(res_path)) {
  stop("Missing DESeq2 result table: ", res_path)
}

dir.create("results/task12", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

res <- read.delim(res_path, stringsAsFactors = FALSE)

required_cols <- c("gene_id", "gene_name", "baseMean", "log2FoldChange", "stat", "pvalue", "padj")
missing_cols <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop("Missing columns in DESeq2 result table: ", paste(missing_cols, collapse = ", "))
}

message("Loaded DESeq2 results: ", nrow(res), " genes")

# Remove Ensembl version number, e.g. ENSG00000162614.21 -> ENSG00000162614
res$ensembl_id <- sub("\\..*$", "", res$gene_id)

# Background universe: all genes tested by DESeq2 and mapped to Entrez
universe_ensembl <- unique(res$ensembl_id)

message("Mapping Ensembl IDs to Entrez IDs...")

id_map <- bitr(
  universe_ensembl,
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

id_map <- unique(id_map)

write.table(
  id_map,
  file = "results/task12/ensembl_to_entrez_mapping.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Mapped genes: ", length(unique(id_map$ENSEMBL)), " Ensembl IDs")
message("Mapped Entrez IDs: ", length(unique(id_map$ENTREZID)))

res_mapped <- merge(
  res,
  id_map,
  by.x = "ensembl_id",
  by.y = "ENSEMBL",
  all.x = FALSE,
  sort = FALSE
)

message("DESeq2 result rows after ID mapping: ", nrow(res_mapped))

# -----------------------------
# 1. Define up/down gene lists
# -----------------------------

sig <- subset(res_mapped, !is.na(padj) & padj < 0.05)
sig_lfc1 <- subset(sig, abs(log2FoldChange) >= 1)

up <- subset(sig_lfc1, log2FoldChange > 0)
down <- subset(sig_lfc1, log2FoldChange < 0)

up_entrez <- unique(up$ENTREZID)
down_entrez <- unique(down$ENTREZID)
universe_entrez <- unique(res_mapped$ENTREZID)

message("Upregulated Entrez genes: ", length(up_entrez))
message("Downregulated Entrez genes: ", length(down_entrez))
message("Universe Entrez genes: ", length(universe_entrez))

write.table(
  up[, c("gene_id", "gene_name", "ENSEMBL_SYMBOL" = "SYMBOL", "ENTREZID", "baseMean", "log2FoldChange", "padj")],
  file = "results/task12/upregulated_genes_for_go.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  down[, c("gene_id", "gene_name", "ENSEMBL_SYMBOL" = "SYMBOL", "ENTREZID", "baseMean", "log2FoldChange", "padj")],
  file = "results/task12/downregulated_genes_for_go.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 2. ORA: GO Biological Process enrichment
# -----------------------------

message("Running GO BP over-representation analysis for upregulated genes...")

ego_up <- enrichGO(
  gene = up_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

message("Running GO BP over-representation analysis for downregulated genes...")

ego_down <- enrichGO(
  gene = down_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_up_df <- as.data.frame(ego_up)
ego_down_df <- as.data.frame(ego_down)

write.table(
  ego_up_df,
  file = "results/task12/go_bp_enrichment_upregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ego_down_df,
  file = "results/task12/go_bp_enrichment_downregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Dotplots
if (nrow(ego_up_df) > 0) {
  pdf("results/figures/task12_go_bp_up_dotplot.pdf", width = 8, height = 6)
  print(dotplot(ego_up, showCategory = 20) + ggtitle("GO BP enrichment: upregulated genes"))
  dev.off()
}

if (nrow(ego_down_df) > 0) {
  pdf("results/figures/task12_go_bp_down_dotplot.pdf", width = 8, height = 6)
  print(dotplot(ego_down, showCategory = 20) + ggtitle("GO BP enrichment: downregulated genes"))
  dev.off()
}

# -----------------------------
# 3. GSEA: ranked gene list
# -----------------------------

message("Preparing ranked gene list for GSEA...")

# Use DESeq2 Wald statistic as ranking metric.
# Remove NA stats and duplicate Entrez IDs.
rank_df <- res_mapped[!is.na(res_mapped$stat) & !is.na(res_mapped$ENTREZID), ]

# If multiple Ensembl IDs map to the same Entrez ID, keep the one with largest absolute statistic.
rank_df <- rank_df %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list <- rank_df$stat
names(gene_list) <- rank_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

write.table(
  data.frame(ENTREZID = names(gene_list), stat = as.numeric(gene_list)),
  file = "results/task12/gsea_ranked_gene_list.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Ranked genes for GSEA: ", length(gene_list))

message("Running GO BP GSEA...")

gsea_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_bp_df <- as.data.frame(gsea_bp)

write.table(
  gsea_bp_df,
  file = "results/task12/gsea_go_bp.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# GSEA dotplot
if (nrow(gsea_bp_df) > 0) {
  pdf("results/figures/task12_gsea_go_bp_dotplot.pdf", width = 8, height = 7)
  print(dotplot(gsea_bp, showCategory = 20, split = ".sign") +
          facet_grid(. ~ .sign) +
          ggtitle("GSEA GO BP: dex treated vs untreated"))
  dev.off()

  # Enrichment plots for top positive and negative terms
  gsea_pos <- subset(gsea_bp_df, NES > 0)
  gsea_neg <- subset(gsea_bp_df, NES < 0)

  if (nrow(gsea_pos) > 0) {
    top_pos_id <- gsea_pos$ID[which.min(gsea_pos$p.adjust)]
    pdf("results/figures/task12_gsea_top_positive_enrichment_plot.pdf", width = 8, height = 6)
    print(gseaplot2(gsea_bp, geneSetID = top_pos_id, title = gsea_bp_df$Description[gsea_bp_df$ID == top_pos_id][1]))
    dev.off()
  }

  if (nrow(gsea_neg) > 0) {
    top_neg_id <- gsea_neg$ID[which.min(gsea_neg$p.adjust)]
    pdf("results/figures/task12_gsea_top_negative_enrichment_plot.pdf", width = 8, height = 6)
    print(gseaplot2(gsea_bp, geneSetID = top_neg_id, title = gsea_bp_df$Description[gsea_bp_df$ID == top_neg_id][1]))
    dev.off()
  }
}

# -----------------------------
# 4. Summary
# -----------------------------

summary_df <- data.frame(
  item = c(
    "genes_in_deseq2_result",
    "mapped_ensembl_ids",
    "mapped_entrez_ids",
    "upregulated_entrez_for_ora",
    "downregulated_entrez_for_ora",
    "go_bp_up_terms",
    "go_bp_down_terms",
    "ranked_genes_for_gsea",
    "gsea_go_bp_terms",
    "gsea_go_bp_positive_terms",
    "gsea_go_bp_negative_terms"
  ),
  value = c(
    nrow(res),
    length(unique(id_map$ENSEMBL)),
    length(unique(id_map$ENTREZID)),
    length(up_entrez),
    length(down_entrez),
    nrow(ego_up_df),
    nrow(ego_down_df),
    length(gene_list),
    nrow(gsea_bp_df),
    sum(gsea_bp_df$NES > 0, na.rm = TRUE),
    sum(gsea_bp_df$NES < 0, na.rm = TRUE)
  )
)

write.table(
  summary_df,
  file = "results/task12/task12_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Task 12 finished.")
message("Main outputs:")
message("  results/task12/go_bp_enrichment_upregulated.tsv")
message("  results/task12/go_bp_enrichment_downregulated.tsv")
message("  results/task12/gsea_go_bp.tsv")
message("  results/figures/task12_go_bp_up_dotplot.pdf")
message("  results/figures/task12_go_bp_down_dotplot.pdf")
message("  results/figures/task12_gsea_go_bp_dotplot.pdf")
message("  results/task12/task12_summary.tsv")
