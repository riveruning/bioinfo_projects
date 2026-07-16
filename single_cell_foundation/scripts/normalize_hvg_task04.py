import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

print("Task 04: normalization, log transform, and highly variable genes")

# -----------------------------
# 1. Paths
# -----------------------------
input_h5ad = Path("data/pbmc3k_qc_filtered.h5ad")
data_dir = Path("data")
results_dir = Path("results/task04")
figures_dir = Path("figures/task04")

results_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

if not input_h5ad.exists():
    raise FileNotFoundError(f"Missing input file: {input_h5ad}")

# -----------------------------
# 2. Load QC-filtered AnnData
# -----------------------------
adata = sc.read_h5ad(input_h5ad)

print("\nInput AnnData:")
print(adata)

n_cells, n_genes = adata.shape
print(f"Cells: {n_cells}")
print(f"Genes: {n_genes}")

# -----------------------------
# 3. Preserve raw counts
# -----------------------------
# adata.X currently contains raw counts after QC filtering.
# Save a copy in layers["counts"] before normalization.
adata.layers["counts"] = adata.X.copy()

raw_total_counts = np.asarray(adata.X.sum(axis=1)).ravel()

# -----------------------------
# 4. Normalize total counts per cell
# -----------------------------
# Normalize every cell to the same total count scale.
sc.pp.normalize_total(adata, target_sum=1e4)

normalized_total_counts = np.asarray(adata.X.sum(axis=1)).ravel()

# Store normalized but not log-transformed values.
adata.layers["normalized"] = adata.X.copy()

# -----------------------------
# 5. Log transform
# -----------------------------
sc.pp.log1p(adata)

# After this step, adata.X contains log-normalized expression values.

# -----------------------------
# 6. Find highly variable genes
# -----------------------------
# This follows the classic Scanpy PBMC3k tutorial-style thresholds.
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

n_hvg = int(adata.var["highly_variable"].sum())
print(f"\nHighly variable genes: {n_hvg}")

# -----------------------------
# 7. Save processed AnnData
# -----------------------------
output_h5ad = data_dir / "pbmc3k_normalized_hvg.h5ad"
adata.write(output_h5ad)
print(f"Saved normalized AnnData to: {output_h5ad}")

# -----------------------------
# 8. Save summary tables
# -----------------------------
summary = pd.DataFrame({
    "item": [
        "cells",
        "genes_after_qc",
        "highly_variable_genes",
        "normalization_target_sum",
        "raw_total_counts_mean",
        "raw_total_counts_median",
        "normalized_total_counts_mean",
        "normalized_total_counts_median",
        "hvg_min_mean",
        "hvg_max_mean",
        "hvg_min_disp"
    ],
    "value": [
        n_cells,
        n_genes,
        n_hvg,
        10000,
        float(np.mean(raw_total_counts)),
        float(np.median(raw_total_counts)),
        float(np.mean(normalized_total_counts)),
        float(np.median(normalized_total_counts)),
        0.0125,
        3,
        0.5
    ]
})

summary.to_csv(results_dir / "task04_normalization_hvg_summary.tsv", sep="\t", index=False)

# Save all gene-level HVG statistics.
hvg_columns = [
    "gene_ids",
    "mt",
    "n_cells",
    "highly_variable",
    "means",
    "dispersions",
    "dispersions_norm"
]

available_hvg_columns = [c for c in hvg_columns if c in adata.var.columns]

adata.var[available_hvg_columns].to_csv(
    results_dir / "pbmc3k_hvg_table.tsv",
    sep="\t"
)

# Save only highly variable genes.
adata.var.loc[adata.var["highly_variable"], available_hvg_columns].to_csv(
    results_dir / "pbmc3k_highly_variable_genes.tsv",
    sep="\t"
)

# Save top 30 HVGs by normalized dispersion.
top_hvg = (
    adata.var.loc[adata.var["highly_variable"], available_hvg_columns]
    .sort_values("dispersions_norm", ascending=False)
    .head(30)
)

top_hvg.to_csv(results_dir / "pbmc3k_top30_hvg.tsv", sep="\t")

# -----------------------------
# 9. Plots
# -----------------------------

# Raw total counts distribution
plt.figure(figsize=(6, 4))
plt.hist(raw_total_counts, bins=50)
plt.xlabel("Raw total counts per cell")
plt.ylabel("Number of cells")
plt.title("Before normalization")
plt.tight_layout()
plt.savefig(figures_dir / "raw_total_counts_hist.png", dpi=150)
plt.close()

# Normalized total counts distribution
plt.figure(figsize=(6, 4))
plt.hist(normalized_total_counts, bins=1, range=(9999, 10001))
plt.xlabel("Normalized total counts per cell")
plt.ylabel("Number of cells")
plt.title("After normalize_total target_sum=10000")
plt.tight_layout()
plt.savefig(figures_dir / "normalized_total_counts_hist.png", dpi=150)
plt.close()

# HVG mean-dispersion plot
plt.figure(figsize=(6, 5))

not_hvg = ~adata.var["highly_variable"]
is_hvg = adata.var["highly_variable"]

plt.scatter(
    adata.var.loc[not_hvg, "means"],
    adata.var.loc[not_hvg, "dispersions_norm"],
    s=5,
    alpha=0.4,
    label="not HVG"
)

plt.scatter(
    adata.var.loc[is_hvg, "means"],
    adata.var.loc[is_hvg, "dispersions_norm"],
    s=5,
    alpha=0.6,
    label="HVG"
)

plt.axhline(0.5, linestyle="--", linewidth=1)
plt.xlabel("Mean expression")
plt.ylabel("Normalized dispersion")
plt.title("Highly variable gene selection")
plt.legend(markerscale=3)
plt.tight_layout()
plt.savefig(figures_dir / "hvg_mean_dispersion_plot.png", dpi=150)
plt.close()

print("\nSummary:")
print(summary)

print("\nTop 30 HVGs:")
print(top_hvg.head(30))

print("\nMain outputs:")
print(f"  {output_h5ad}")
print(f"  {results_dir / 'task04_normalization_hvg_summary.tsv'}")
print(f"  {results_dir / 'pbmc3k_hvg_table.tsv'}")
print(f"  {results_dir / 'pbmc3k_highly_variable_genes.tsv'}")
print(f"  {results_dir / 'pbmc3k_top30_hvg.tsv'}")
print(f"  {figures_dir}")
print("Task 04 finished.")
