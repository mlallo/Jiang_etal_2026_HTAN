import os
from pathlib import Path
import scanpy as sc
import hotspot
from scipy import sparse

## Input The Main Directory ##
main_directory = Path("/path/to/Jiang_etal_2026_HTAN")

idir = main_directory / "output" / "processed_adata"
hs_dir = main_directory / "output" / "figures" / "figures_extended_data" / "figure_1g"
hs_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(idir / "processed_adata.h5ad")

adata = adata[~adata.obs['Cell Type'].isin(['Primary Tumor', 'Metastasis'])].copy()

sc.pp.filter_genes(adata, min_cells=1)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=False,
    flavor="seurat_v3",
    layer='raw'
)

adata.var.loc['ZFP36L2', 'highly_variable'] = True

bad_genes = adata.var_names.str.contains(
    r"^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT"
)
adata.var.loc[bad_genes, 'highly_variable'] = False

# Subset the AnnData object to highly variable genes
adata = adata[:, adata.var.highly_variable].copy()

print("Remaining Cell Types after filtering:", adata.obs['Cell Type'].unique().tolist())

sc.pp.neighbors(adata, n_neighbors=50, n_pcs=adata.obsm['X_pca'].shape[1])

adata.layers['counts'] = adata.raw[adata.obs_names, adata.var_names].X

counts = adata.layers['counts']
if not sparse.issparse(counts):
    counts = sparse.csr_matrix(counts)

nonzero_gene_mask = counts.sum(axis=0).A1 > 0
adata = adata[:, nonzero_gene_mask].copy()

print(f"Number of genes with non-zero counts: {adata.n_vars}")

print(f"Number of genes for hotspot: {adata.n_vars}")

# Initialize Hotspot
hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    distances_obsp_key="distances",
    umi_counts_obs_key="total_counts"
)

# Create KNN graph
hs.create_knn_graph(weighted_graph=False, n_neighbors=15)

# Compute autocorrelations
hs_results = hs.compute_autocorrelations(jobs=30)
#hs_results.to_csv(os.path.join(hs_dir, 'hotspot.informative_genes.tsv'), sep='\t')

# Select significant genes
hs_genes = hs_results.loc[hs_results['FDR'] < 0.05].index

# Compute local correlations
local_correlations = hs.compute_local_correlations(hs_genes, jobs=30)
local_correlations.to_csv(os.path.join(hs_dir, 'hotspot.autocorrelations_normal.tsv'), sep='\t')
