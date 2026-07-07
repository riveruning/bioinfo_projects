# Inspect Salmon quantification result for SRR1039508

quant <- read.delim("results/salmon/SRR1039508/quant.sf")

cat("Number of transcripts:", nrow(quant), "\n")
cat("Total estimated NumReads:", sum(quant$NumReads), "\n")
cat("Detected transcripts with TPM > 0:", sum(quant$TPM > 0), "\n")

top20 <- quant[order(-quant$TPM), ][1:20, ]

cat("\nTop 20 transcripts by TPM:\n")
print(top20[, c("Name", "Length", "EffectiveLength", "TPM", "NumReads")])
