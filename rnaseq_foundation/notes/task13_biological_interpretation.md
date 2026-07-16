# Task 13: Biological Interpretation of DESeq2 and Enrichment Results

## 1. Task Goal

The goal of this task was to summarize the main biological signals from the DESeq2 differential expression results and GO/GSEA enrichment analysis.

This task focused on identifying:

- key dexamethasone-responsive genes
- major biological processes affected by treatment
- candidate genes for potential experimental validation

## 2. Input Files

Main input files:

- results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv
- results/task11/candidate_genes_padj_0.05_abs_log2fc_1_basemean_100.tsv
- results/task12/go_bp_enrichment_upregulated.tsv
- results/task12/go_bp_enrichment_downregulated.tsv
- results/task12/gsea_go_bp.tsv

## 3. Main Output Files

Task 13 generated:

- results/task13/key_gene_interpretation_table.tsv
- results/task13/validation_candidate_genes.tsv
- results/task13/pathway_interpretation_summary.tsv
- results/task13/gsea_terms_suitable_for_report.tsv
- results/task13/task13_summary.tsv

Summary:

- key genes found: 12
- validation candidates: 8
- top GO terms summarized: 30
- GSEA report terms: 40

## 4. Key Genes

Important upregulated genes included:

- ZBTB16
- DUSP1
- PER1
- CRISPLD2
- GPX3
- ADAMTS1
- SERPINA3
- SAMHD1
- MT2A
- SORT1
- NEXN

Important downregulated genes included:

- VCAM1

## 5. Candidate Genes for Validation

The main validation candidates were:

- ZBTB16: strong glucocorticoid-responsive upregulated gene
- DUSP1: stress/MAPK-related glucocorticoid-responsive gene
- PER1: hormone/circadian-related response gene
- CRISPLD2: known airway glucocorticoid-responsive candidate from the original dataset paper
- VCAM1: strongly downregulated adhesion/inflammatory marker
- ADAMTS1: extracellular matrix remodeling candidate
- GPX3: oxidative stress or antioxidant-response candidate
- SERPINA3: secreted/inflammatory-response candidate

## 6. Approximate Fold Changes

Selected genes showed strong treatment effects:

- ZBTB16: log2FC about 6.10, approximately 68-fold upregulated
- DUSP1: log2FC about 2.97, approximately 7.8-fold upregulated
- PER1: log2FC about 3.03, approximately 8.2-fold upregulated
- CRISPLD2: log2FC about 2.61, approximately 6.1-fold upregulated
- VCAM1: log2FC about -3.59, approximately 12-fold downregulated

For downregulated genes, fold-change values should be described as lower expression rather than negative expression.

## 7. Main Upregulated Biological Processes

Upregulated genes were enriched for:

- muscle system process
- cell-substrate adhesion
- muscle adaptation
- angiogenesis
- cellular response to glucocorticoid stimulus
- cellular response to corticosteroid stimulus
- reactive oxygen species metabolic process
- extracellular matrix organization
- extracellular structure organization
- regulation of blood pressure

These suggest that dexamethasone treatment affects glucocorticoid response, muscle-related programs, extracellular matrix organization, adhesion, and oxidative stress-related processes.

## 8. Main Downregulated Biological Processes

Downregulated genes were enriched for GO terms related to:

- regulation of nervous system development
- axon development
- axon guidance
- animal organ morphogenesis
- chemotaxis
- cell migration
- cell motility
- MAPK cascade regulation

These terms should not be interpreted literally as neuronal differentiation in airway smooth muscle cells.

A more cautious interpretation is that dexamethasone treatment downregulates genes involved in developmental signaling, cell guidance, morphogenesis, migration, and chemotaxis-related processes.

## 9. GSEA Interpretation

GSEA supported activated processes related to:

- cellular response to hormone stimulus
- glucose metabolism
- D-glucose import
- carbohydrate metabolism
- energy metabolism
- cell-substrate adhesion
- muscle system process
- response to oxidative stress

Some suppressed GSEA terms included developmental or morphogenesis-related processes.

GO BP GSEA was useful as a global trend check, but the main interpretation should rely more heavily on ORA results and biologically interpretable terms.

## 10. Overall Biological Interpretation

Dexamethasone treatment caused broad transcriptional remodeling in airway smooth muscle cells.

The treatment induced glucocorticoid-responsive genes such as ZBTB16, DUSP1, PER1, and CRISPLD2.

It also affected extracellular matrix organization, cell-substrate adhesion, muscle-related processes, oxidative stress response, and glucose/energy metabolism.

At the same time, dexamethasone downregulated VCAM1 and gene programs associated with migration, chemotaxis, morphogenesis, and developmental signaling.

Overall, these results are consistent with a strong glucocorticoid response and suggest that dexamethasone alters inflammatory/adhesion-related and tissue-remodeling-related programs in airway smooth muscle cells.

## 11. Next Step

The next step is Task 14: write the final project report and README.

The final report should summarize:

- dataset background
- analysis workflow
- QC results
- DESeq2 results
- key figures
- enrichment analysis
- biological interpretation
- limitations
- possible experimental validation
