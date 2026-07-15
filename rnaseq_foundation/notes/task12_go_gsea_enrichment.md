# Task 12: GO Enrichment and GSEA

## 1. Task Goal

The goal of this task was to interpret the DESeq2 differential expression results at the pathway and biological process level.

Instead of only looking at individual differentially expressed genes, this task asked:

- Which biological processes are enriched among upregulated genes?
- Which biological processes are enriched among downregulated genes?
- Which biological processes are globally activated or suppressed across the full ranked gene list?

## 2. Input Files

The main input file was:

- results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv

This file contains DESeq2 differential expression results for dexamethasone-treated samples versus untreated samples.

Important columns used:

- gene_id
- gene_name
- baseMean
- log2FoldChange
- stat
- pvalue
- padj

## 3. ID Mapping

The DESeq2 result table used Ensembl gene IDs with version numbers, such as:

ENSG00000162614.21

For GO enrichment analysis, Ensembl version numbers were removed:

ENSG00000162614.21 -> ENSG00000162614

Then Ensembl IDs were mapped to Entrez IDs using org.Hs.eg.db.

Summary:

- genes in DESeq2 result: 20407
- mapped Ensembl IDs: 16902
- mapped Entrez IDs: 16571

The mapped Entrez IDs were used for GO enrichment and GSEA.

## 4. ORA: GO Biological Process Enrichment

Over-representation analysis, ORA, was performed separately for upregulated and downregulated genes.

The input gene lists were selected using:

- padj < 0.05
- absolute log2FoldChange >= 1

Summary:

- upregulated Entrez genes used for ORA: 485
- downregulated Entrez genes used for ORA: 433
- enriched GO BP terms for upregulated genes: 77
- enriched GO BP terms for downregulated genes: 418

Output files:

- results/task12/go_bp_enrichment_upregulated.tsv
- results/task12/go_bp_enrichment_downregulated.tsv
- results/task12/upregulated_genes_for_go.tsv
- results/task12/downregulated_genes_for_go.tsv

Figures:

- results/figures/task12_go_bp_up_dotplot.pdf
- results/figures/task12_go_bp_down_dotplot.pdf

## 5. Upregulated Gene Enrichment

Upregulated genes were enriched for biological processes related to:

- muscle system process
- cell-substrate adhesion
- regulation of system process
- muscle adaptation
- angiogenesis
- cellular response to glucocorticoid stimulus
- regulation of blood pressure
- reactive oxygen species metabolic process
- extracellular matrix organization
- extracellular structure organization
- response to hypoxia
- cellular response to corticosteroid stimulus

These results are broadly consistent with dexamethasone-induced transcriptional remodeling in airway smooth muscle cells.

Important terms include:

- cellular response to glucocorticoid stimulus
- cellular response to corticosteroid stimulus

These terms directly support the presence of a glucocorticoid-response signal in the dataset.

## 6. Downregulated Gene Enrichment

Downregulated genes were enriched for biological processes related to:

- positive regulation of nervous system development
- regulation of nervous system development
- axon development
- axonogenesis
- axon guidance
- chemotaxis
- taxis
- positive regulation of cell migration
- positive regulation of cell motility
- circulatory system process
- blood circulation
- positive regulation of MAPK cascade

Some GO terms mention nervous system development or axon guidance.

These should not be interpreted literally as airway smooth muscle cells becoming neuronal cells.

Many genes annotated to nervous system development or axon guidance, such as signaling, migration, guidance, and morphogenesis-related genes, can also participate in broader processes such as:

- cell migration
- chemotaxis
- morphogenesis
- cell-cell signaling
- tissue remodeling

A more cautious interpretation is:

Downregulated genes were enriched for developmental signaling, chemotaxis, cell migration, and morphogenesis-related biological processes.

## 7. GSEA: GO Biological Process

GSEA was performed using the full ranked gene list.

The ranking metric was the DESeq2 Wald statistic.

Interpretation:

- NES > 0: the GO term is globally enriched toward genes upregulated in dexamethasone-treated samples
- NES < 0: the GO term is globally enriched toward genes downregulated in dexamethasone-treated samples

Summary:

- ranked genes used for GSEA: 16571
- significant GSEA GO BP terms: 229
- positive GSEA terms: 185
- negative GSEA terms: 44

Output files:

- results/task12/gsea_ranked_gene_list.tsv
- results/task12/gsea_go_bp.tsv

Figures:

- results/figures/task12_gsea_go_bp_dotplot.pdf
- results/figures/task12_gsea_top_positive_enrichment_plot.pdf
- results/figures/task12_gsea_top_negative_enrichment_plot.pdf

## 8. GSEA Activated Processes

GSEA activated terms included processes related to:

- cellular response to peptide hormone stimulus
- response to oxidative stress
- generation of precursor metabolites and energy
- energy derivation by oxidation of organic compounds
- carbohydrate metabolic process
- glucose metabolic process
- D-glucose import
- hexose transmembrane transport
- cell-substrate adhesion
- muscle system process

These results suggest that dexamethasone treatment is associated with broad changes in hormone-response-related processes, energy metabolism, glucose transport, oxidative stress response, and adhesion/muscle-related processes.

The top positive enrichment plot showed that genes in the term:

- cellular response to peptide hormone stimulus

were enriched toward the dexamethasone-treated upregulated side of the ranked gene list.

However, this term should be interpreted cautiously. Dexamethasone is not a peptide hormone. The term is broader and should be described as hormone-response-related rather than as a direct peptide hormone response.

## 9. GSEA Suppressed Processes

GSEA suppressed terms included processes related to:

- hair follicle morphogenesis
- epidermis morphogenesis
- sensory organ morphogenesis
- regulation of nervous system development
- forebrain development
- cell fate specification
- cell fate commitment

The top negative enrichment plot showed:

- hair follicle morphogenesis

as the strongest suppressed term.

This should not be interpreted literally as dexamethasone suppressing hair follicle formation in airway smooth muscle cells.

A more cautious interpretation is:

Some developmental, morphogenesis, and cell fate-related gene sets were globally shifted toward the untreated or downregulated side of the ranked list after dexamethasone treatment.

## 10. Initial Biological Interpretation

Overall, Task 12 suggests that dexamethasone treatment causes broad transcriptional remodeling in airway smooth muscle cells.

The upregulated direction is associated with:

- glucocorticoid/corticosteroid response
- muscle-related processes
- extracellular matrix organization
- cell-substrate adhesion
- angiogenesis-related processes
- oxidative stress response
- hypoxia response
- glucose and energy metabolism

The downregulated direction is associated with:

- chemotaxis
- cell migration
- cell motility
- morphogenesis
- developmental signaling
- nervous-system-development-like GO annotations

The most important interpretation is not that every GO term should be read literally.

Instead, GO enrichment provides a functional overview suggesting that dexamethasone affects hormone response, metabolism, extracellular matrix/adhesion, and migration/developmental signaling programs.

## 11. Important Cautions

GO terms are annotations, not direct experimental proof.

Some GO terms can sound tissue-specific or development-specific, such as:

- axon guidance
- nervous system development
- hair follicle morphogenesis
- forebrain development

These terms may appear because the same genes are involved in broader biological processes such as cell migration, guidance, morphogenesis, and signaling.

Therefore, these results should be interpreted as functional hints, not final mechanistic conclusions.

Further interpretation should combine:

- top differentially expressed genes
- GO/GSEA results
- known dexamethasone biology
- airway smooth muscle cell context
- literature evidence
- possible experimental validation

## 12. Next Step

The next step is Task 13: biological interpretation of key genes and pathways.

This will include:

- interpreting major enriched biological processes
- checking known dexamethasone-responsive genes
- summarizing candidate genes such as ZBTB16, DUSP1, PER1, VCAM1, CRISPLD2, GPX3, ADAMTS1, and SERPINA3
- selecting possible genes for qPCR or functional validation
