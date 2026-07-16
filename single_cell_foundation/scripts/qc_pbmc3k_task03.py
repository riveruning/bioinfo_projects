import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

print("Task 03: PBMC3k quality control")

# -----------------------------
# 1. Paths
# -----------------------------
input_h5ad = Path("data/pbmc3k_raw.h5ad")
results_dir = Path("results/task03")
figures_dir = Path("figures/task03")
data_dir = Path("data")

results_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

if not input_h5ad.exists():
    raise FileNotFoundError(f"Missing input file: {input_h5ad}")

# -----------------------------
# 2. Load raw AnnData
# -----------------------------
adata = sc.read_h5ad(input_h5ad)

print("\nRaw AnnData:")
print(adata)

n_cells_raw, n_genes_raw = adata.shape
print(f"Raw cells: {n_cells_raw}")
print(f"Raw genes: {n_genes_raw}")

# -----------------------------
# 3. Mark mitochondrial genes
# -----------------------------
# Human mitochondrial genes usually start with "MT-".
adata.var["mt"] = adata.var_names.str.startswith("MT-")

n_mt_genes = int(adata.var["mt"].sum())
print(f"Mitochondrial genes detected by MT- prefix: {n_mt_genes}")

# -----------------------------
# 4. Calculate QC metrics
# -----------------------------
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("\nobs columns after QC:")
print(list(adata.obs.columns))

# -----------------------------
# 5. Save raw QC metrics
# -----------------------------
qc_cols = [
    "n_genes_by_counts",
    "total_counts",
    "total_counts_mt",
    "pct_counts_mt"
]

adata.obs[qc_cols].to_csv(
    results_dir / "pbmc3k_cell_qc_metrics_raw.tsv",
    sep="\t"
)

qc_summary_raw = adata.obs[qc_cols].describe().T
qc_summary_raw.to_csv(
    results_dir / "pbmc3k_qc_summary_raw.tsv",
    sep="\t"
)

# -----------------------------
# 6. Plot QC distributions before filtering
# -----------------------------
for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
    plt.figure(figsize=(6, 4))
    plt.hist(adata.obs[col], bins=50)
    plt.xlabel(col)
    plt.ylabel("Number of cells")
    plt.title(f"PBMC3k raw QC: {col}")
    plt.tight_layout()
    plt.savefig(figures_dir / f"raw_{col}_hist.png", dpi=150)
    plt.close()

plt.figure(figsize=(5, 4))
plt.scatter(
    adata.obs["total_counts"],
    adata.obs["pct_counts_mt"],
    s=6,
    alpha=0.6
)
plt.xlabel("total_counts")
plt.ylabel("pct_counts_mt")
plt.title("Raw cells: total counts vs mitochondrial percentage")
plt.tight_layout()
plt.savefig(figures_dir / "raw_total_counts_vs_pct_counts_mt.png", dpi=150)
plt.close()

plt.figure(figsize=(5, 4))
plt.scatter(
    adata.obs["total_counts"],
    adata.obs["n_genes_by_counts"],
    s=6,
    alpha=0.6
)
plt.xlabel("total_counts")
plt.ylabel("n_genes_by_counts")
plt.title("Raw cells: total counts vs detected genes")
plt.tight_layout()
plt.savefig(figures_dir / "raw_total_counts_vs_n_genes.png", dpi=150)
plt.close()

# -----------------------------
# 7. Filter cells and genes
# -----------------------------
# Classic PBMC3k filtering logic:
# - keep cells with at least 200 detected genes
# - keep genes detected in at least 3 cells
# - remove cells with too many detected genes, possible doublets
# - remove cells with high mitochondrial percentage

adata_filtered = adata.copy()

sc.pp.filter_cells(adata_filtered, min_genes=200)
sc.pp.filter_genes(adata_filtered, min_cells=3)

adata_filtered = adata_filtered[
    adata_filtered.obs["n_genes_by_counts"] < 2500,
    :
].copy()

adata_filtered = adata_filtered[
    adata_filtered.obs["pct_counts_mt"] < 5,
    :
].copy()

n_cells_filtered, n_genes_filtered = adata_filtered.shape

print("\nFiltered AnnData:")
print(adata_filtered)
print(f"Filtered cells: {n_cells_filtered}")
print(f"Filtered genes: {n_genes_filtered}")

# -----------------------------
# 8. Save filtered QC metrics and object
# -----------------------------
adata_filtered.obs[qc_cols].to_csv(
    results_dir / "pbmc3k_cell_qc_metrics_filtered.tsv",
    sep="\t"
)

qc_summary_filtered = adata_filtered.obs[qc_cols].describe().T
qc_summary_filtered.to_csv(
    results_dir / "pbmc3k_qc_summary_filtered.tsv",
    sep="\t"
)

adata_filtered.write(data_dir / "pbmc3k_qc_filtered.h5ad")

# -----------------------------
# 9. Plot QC distributions after filtering
# -----------------------------
for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
    plt.figure(figsize=(6, 4))
    plt.hist(adata_filtered.obs[col], bins=50)
    plt.xlabel(col)
    plt.ylabel("Number of cells")
    plt.title(f"PBMC3k filtered QC: {col}")
    plt.tight_layout()
    plt.savefig(figures_dir / f"filtered_{col}_hist.png", dpi=150)
    plt.close()

# -----------------------------
# 10. Save filtering summary
# -----------------------------
summary = pd.DataFrame({
    "item": [
        "raw_cells",
        "raw_genes",
        "filtered_cells",
        "filtered_genes",
        "removed_cells",
        "removed_genes",
        "mitochondrial_genes_detected",
        "cell_filter_min_genes",
        "cell_filter_max_genes",
        "cell_filter_max_pct_counts_mt",
        "gene_filter_min_cells"
    ],
    "value": [
        n_cells_raw,
        n_genes_raw,
        n_cells_filtered,
        n_genes_filtered,
        n_cells_raw - n_cells_filtered,
        n_genes_raw - n_genes_filtered,
        n_mt_genes,
        200,
        2500,
        5,
        3
    ]
})

summary.to_csv(results_dir / "pbmc3k_qc_filtering_summary.tsv", sep="\t", index=False)

print("\nFiltering summary:")
print(summary)

print("\nMain outputs:")
print(f"  {results_dir / 'pbmc3k_qc_filtering_summary.tsv'}")
print(f"  {results_dir / 'pbmc3k_qc_summary_raw.tsv'}")
print(f"  {results_dir / 'pbmc3k_qc_summary_filtered.tsv'}")
print(f"  {data_dir / 'pbmc3k_qc_filtered.h5ad'}")
print(f"  {figures_dir}")
print("Task 03 finished.")
