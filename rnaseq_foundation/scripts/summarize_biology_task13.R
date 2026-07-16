suppressPackageStartupMessages({
  library(dplyr)
})

message("Task 13: biological interpretation summary")

res_path <- "results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv"
candidate_path <- "results/task11/candidate_genes_padj_0.05_abs_log2fc_1_basemean_100.tsv"
go_up_path <- "results/task12/go_bp_enrichment_upregulated.tsv"
go_down_path <- "results/task12/go_bp_enrichment_downregulated.tsv"
gsea_path <- "results/task12/gsea_go_bp.tsv"

if (!file.exists(res_path)) stop("Missing: ", res_path)
if (!file.exists(candidate_path)) stop("Missing: ", candidate_path)
if (!file.exists(go_up_path)) stop("Missing: ", go_up_path)
if (!file.exists(go_down_path)) stop("Missing: ", go_down_path)
if (!file.exists(gsea_path)) stop("Missing: ", gsea_path)

dir.create("results/task13", recursive = TRUE, showWarnings = FALSE)

res <- read.delim(res_path, stringsAsFactors = FALSE)
candidate <- read.delim(candidate_path, stringsAsFactors = FALSE)
go_up <- read.delim(go_up_path, stringsAsFactors = FALSE)
go_down <- read.delim(go_down_path, stringsAsFactors = FALSE)
gsea <- read.delim(gsea_path, stringsAsFactors = FALSE)

# -----------------------------
# 1. Key genes for interpretation
# -----------------------------

key_genes <- c(
  "ZBTB16", "DUSP1", "PER1", "VCAM1",
  "CRISPLD2", "GPX3", "ADAMTS1", "SERPINA3",
  "SAMHD1", "MT2A", "SORT1", "NEXN"
)

key_res <- res %>%
  filter(gene_name %in% key_genes) %>%
  mutate(
    direction = case_when(
      log2FoldChange > 0 ~ "upregulated",
      log2FoldChange < 0 ~ "downregulated",
      TRUE ~ "unchanged"
    ),
    abs_log2FoldChange = abs(log2FoldChange),
    fold_change_approx = ifelse(
      log2FoldChange >= 0,
      2 ^ log2FoldChange,
      -1 * (2 ^ abs(log2FoldChange))
    )
  ) %>%
  arrange(padj)

write.table(
  key_res,
  file = "results/task13/key_gene_interpretation_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 2. Validation candidate genes
# -----------------------------

validation_candidates <- candidate %>%
  filter(
    gene_name %in% c(
      "ZBTB16", "DUSP1", "PER1", "VCAM1",
      "CRISPLD2", "GPX3", "ADAMTS1", "SERPINA3"
    )
  ) %>%
  mutate(
    suggested_validation = case_when(
      gene_name %in% c("ZBTB16", "DUSP1", "PER1", "CRISPLD2") ~ "qPCR: glucocorticoid-responsive upregulated gene",
      gene_name %in% c("VCAM1") ~ "qPCR/flow/IF: downregulated adhesion/inflammatory marker",
      gene_name %in% c("GPX3", "SERPINA3") ~ "qPCR: stress/secreted-response candidate",
      gene_name %in% c("ADAMTS1") ~ "qPCR: extracellular matrix remodeling candidate",
      TRUE ~ "qPCR candidate"
    )
  ) %>%
  arrange(padj)

write.table(
  validation_candidates,
  file = "results/task13/validation_candidate_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 3. Pathway interpretation summary
# -----------------------------

top_go_up <- go_up %>%
  select(ID, Description, GeneRatio, p.adjust, Count) %>%
  head(15) %>%
  mutate(direction = "upregulated genes")

top_go_down <- go_down %>%
  select(ID, Description, GeneRatio, p.adjust, Count) %>%
  head(15) %>%
  mutate(direction = "downregulated genes")

pathway_summary <- bind_rows(top_go_up, top_go_down) %>%
  select(direction, ID, Description, GeneRatio, Count, p.adjust)

write.table(
  pathway_summary,
  file = "results/task13/pathway_interpretation_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 4. GSEA terms suitable for report
# -----------------------------

gsea_report_terms <- gsea %>%
  filter(
    grepl("hormone|oxidative|glucose|carbohydrate|energy|adhesion|muscle|migration|chemotaxis|morphogenesis|MAPK",
          Description,
          ignore.case = TRUE)
  ) %>%
  arrange(p.adjust) %>%
  select(ID, Description, setSize, NES, pvalue, p.adjust, qvalue) %>%
  head(40)

write.table(
  gsea_report_terms,
  file = "results/task13/gsea_terms_suitable_for_report.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 5. Short summary
# -----------------------------

summary_df <- data.frame(
  item = c(
    "key_genes_found",
    "validation_candidates",
    "top_go_terms_summarized",
    "gsea_report_terms"
  ),
  value = c(
    nrow(key_res),
    nrow(validation_candidates),
    nrow(pathway_summary),
    nrow(gsea_report_terms)
  )
)

write.table(
  summary_df,
  file = "results/task13/task13_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Task 13 finished.")
message("Main outputs:")
message("  results/task13/key_gene_interpretation_table.tsv")
message("  results/task13/validation_candidate_genes.tsv")
message("  results/task13/pathway_interpretation_summary.tsv")
message("  results/task13/gsea_terms_suitable_for_report.tsv")
message("  results/task13/task13_summary.tsv")
