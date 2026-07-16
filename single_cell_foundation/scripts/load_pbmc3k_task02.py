import scanpy as sc
from pathlib import Path

print("Task 02: Load PBMC3k dataset")

data_dir = Path("data")
results_dir = Path("results/task02")
data_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

# Load PBMC3k dataset.
# The first run may download the dataset automatically.
adata = sc.datasets.pbmc3k()

print("\nAnnData object:")
print(adata)

n_cells, n_genes = adata.shape
print(f"\nNumber of cells: {n_cells}")
print(f"Number of genes: {n_genes}")

print("\nadata.X:")
print(type(adata.X))
print(adata.X.shape)

print("\nadata.obs head:")
print(adata.obs.head())

print("\nadata.var head:")
print(adata.var.head())

print("\nobs columns:")
print(list(adata.obs.columns))

print("\nvar columns:")
print(list(adata.var.columns))

out_h5ad = data_dir / "pbmc3k_raw.h5ad"
adata.write(out_h5ad)
print(f"\nSaved raw AnnData object to: {out_h5ad}")

summary_path = results_dir / "pbmc3k_basic_summary.tsv"
with open(summary_path, "w") as f:
    f.write("item\tvalue\n")
    f.write(f"n_cells\t{n_cells}\n")
    f.write(f"n_genes\t{n_genes}\n")
    f.write(f"adata_X_type\t{type(adata.X)}\n")
    f.write(f"obs_columns\t{','.join(adata.obs.columns)}\n")
    f.write(f"var_columns\t{','.join(adata.var.columns)}\n")

adata.obs.head(20).to_csv(results_dir / "pbmc3k_obs_head.tsv", sep="\t")
adata.var.head(20).to_csv(results_dir / "pbmc3k_var_head.tsv", sep="\t")

print(f"Saved summary to: {summary_path}")
print("Task 02 finished.")
