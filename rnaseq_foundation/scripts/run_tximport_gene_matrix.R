suppressPackageStartupMessages({
  library(tximport)
})

message("Task 09: tximport gene-level matrix")

samples_path <- "metadata/samples.csv"
tx2gene_full_path <- "metadata/tx2gene_full.tsv"

if (!file.exists(samples_path)) {
  stop("Missing metadata/samples.csv")
}

if (!file.exists(tx2gene_full_path)) {
  stop("Missing metadata/tx2gene_full.tsv")
}

samples <- read.csv(samples_path, stringsAsFactors = FALSE)

if (!"sample_id" %in% colnames(samples)) {
  stop("metadata/samples.csv must contain a sample_id column")
}

sample_ids <- samples$sample_id

files <- file.path("results", "salmon", sample_ids, "quant.sf")
names(files) <- sample_ids

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

message("Found quant.sf files:")
print(files)

tx2gene_full <- read.delim(tx2gene_full_path, stringsAsFactors = FALSE)

required_cols <- c("transcript_id", "gene_id", "gene_name", "gene_type")
missing_cols <- setdiff(required_cols, colnames(tx2gene_full))
if (length(missing_cols) > 0) {
  stop("metadata/tx2gene_full.tsv is missing columns: ", paste(missing_cols, collapse = ", "))
}

# Read one quant.sf file to inspect transcript names used by Salmon
q0 <- read.delim(files[1], stringsAsFactors = FALSE)

if (!"Name" %in% colnames(q0)) {
  stop("quant.sf does not contain a Name column")
}

quant_names <- q0$Name

message("Number of transcript entries in quant.sf: ", length(quant_names))
message("Number of transcript entries in tx2gene_full.tsv: ", nrow(tx2gene_full))

# Case 1: quant.sf Name exactly matches tx2gene transcript_id
exact_overlap <- sum(quant_names %in% tx2gene_full$transcript_id)
exact_rate <- exact_overlap / length(quant_names)

message("Exact transcript ID overlap rate: ", round(exact_rate * 100, 2), "%")

if (exact_rate > 0.95) {
  message("Using exact transcript_id matching.")
  tx2gene_for_salmon <- unique(tx2gene_full[, c("transcript_id", "gene_id")])
} else {
  message("Exact matching is low. Trying to parse GENCODE-style Salmon names.")

  # GENCODE FASTA headers may look like:
  # ENST...|ENSG...|...|GENENAME...
  parsed_transcript_id <- sub("\\|.*$", "", quant_names)

  helper <- data.frame(
    salmon_name = quant_names,
    transcript_id = parsed_transcript_id,
    stringsAsFactors = FALSE
  )

  # First try exact transcript_id after splitting by |
  idx <- match(helper$transcript_id, tx2gene_full$transcript_id)
  gene_id <- tx2gene_full$gene_id[idx]

  # If still missing, try version-stripped transcript IDs
  missing <- is.na(gene_id)
  if (any(missing)) {
    message("Some transcript IDs still unmatched. Trying version-stripped matching.")

    strip_version <- function(x) sub("\\.[0-9]+$", "", x)

    tx2gene_full$transcript_id_no_version <- strip_version(tx2gene_full$transcript_id)
    helper$transcript_id_no_version <- strip_version(helper$transcript_id)

    idx2 <- match(helper$transcript_id_no_version[missing], tx2gene_full$transcript_id_no_version)
    gene_id[missing] <- tx2gene_full$gene_id[idx2]
  }

  matched_rate <- sum(!is.na(gene_id)) / length(gene_id)
  message("Parsed transcript ID matched rate: ", round(matched_rate * 100, 2), "%")

  if (matched_rate < 0.95) {
    stop("Transcript-to-gene matching rate is too low. Stop before tximport.")
  }

  tx2gene_for_salmon <- data.frame(
    transcript_id = helper$salmon_name,
    gene_id = gene_id,
    stringsAsFactors = FALSE
  )

  tx2gene_for_salmon <- tx2gene_for_salmon[!is.na(tx2gene_for_salmon$gene_id), ]
  tx2gene_for_salmon <- unique(tx2gene_for_salmon)
}

dir.create("metadata", showWarnings = FALSE, recursive = TRUE)
write.table(
  tx2gene_for_salmon,
  file = "metadata/tx2gene_for_salmon.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote metadata/tx2gene_for_salmon.tsv")

dir.create("results/tximport", showWarnings = FALSE, recursive = TRUE)

message("Running tximport...")

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene_for_salmon,
  countsFromAbundance = "no",
  dropInfReps = TRUE
)

saveRDS(txi, file = "results/tximport/txi_salmon_gene_level.rds")

write_matrix <- function(mat, path) {
  df <- data.frame(
    gene_id = rownames(mat),
    mat,
    check.names = FALSE
  )

  write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_matrix(txi$counts, "results/tximport/gene_counts.tsv")
write_matrix(txi$abundance, "results/tximport/gene_tpm.tsv")
write_matrix(txi$length, "results/tximport/gene_effective_length.tsv")

gene_annotation <- unique(tx2gene_full[, c("gene_id", "gene_name", "gene_type")])
write.table(
  gene_annotation,
  file = "results/tximport/gene_annotation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- data.frame(
  item = c(
    "samples",
    "genes_in_count_matrix",
    "transcripts_in_quant_sf",
    "transcripts_in_tx2gene_for_salmon"
  ),
  value = c(
    length(sample_ids),
    nrow(txi$counts),
    length(quant_names),
    nrow(tx2gene_for_salmon)
  )
)

write.table(
  summary_df,
  file = "results/tximport/tximport_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Task 09 finished.")
message("Output files:")
message("  results/tximport/txi_salmon_gene_level.rds")
message("  results/tximport/gene_counts.tsv")
message("  results/tximport/gene_tpm.tsv")
message("  results/tximport/gene_effective_length.tsv")
message("  results/tximport/gene_annotation.tsv")
message("  results/tximport/tximport_summary.tsv")
