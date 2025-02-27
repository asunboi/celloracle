import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import anndata as ad
import seaborn as sns

# Initialize the parser
parser = argparse.ArgumentParser(description="Process DEG for GOI, sample, and region.")

# Define the arguments
parser.add_argument('goi', type=str, help='Gene of interest (GOI)')
parser.add_argument('sample', type=int, help='Sample number')
parser.add_argument('region', type=str, help='Region (e.g., ctx)')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
goi = args.goi
sample = args.sample
region = args.region

adata = ad.read_h5ad("/gpfs/home/asun/jinlab/ctxobj.h5ad")

control = ['NonTarget1', 'NonTarget2', 'SafeTarget', 'GFP']
ctrl_adata = adata[adata.obs['assignDS'].isin(control)].copy()

goi_targets = [f"{goi}_1", f"{goi}_2", f"{goi}_3", 'NonTarget1', 'NonTarget2', 'SafeTarget', 'GFP']
goi_adata = adata[adata.obs['assignDS'].isin(goi_targets)].copy()

adata_list = [adata, ctrl_adata, goi_adata]
adata_names = ["peturbseq", "ctrl", f"{goi}"]

for i in range(len(adata_list)):
	current_adata = adata_list[i]
	name = adata_names[i]
	# Only consider genes with more than 1 count
	sc.pp.filter_genes(current_adata, min_counts=1)

	# Normalize gene expression matrix with total UMI count per cell
	sc.pp.normalize_per_cell(current_adata, key_n_counts='nCount_RNA')

	# Get Sean's DEG list
	df = pd.read_csv(f"/gpfs/home/asun/jinlab/tsv/edgeR_LRT_with_sva.{goi}_{sample}.{region}.tsv", sep='\t')
	df_sorted = df.iloc[df['logFC'].abs().argsort()]
	df_sorted.set_index("Unnamed: 0")
	rows_to_add = df_sorted.loc[df_sorted["Unnamed: 0"] == goi]
	DEG = df_sorted.tail(3000)
	DEG = pd.concat([DEG, rows_to_add], ignore_index=True)
	DEG.set_index("Unnamed: 0")
	deg_list = DEG['Unnamed: 0'].unique().tolist()
    
	current_adata = current_adata[:, deg_list]

	# Renormalize after filtering
	sc.pp.normalize_per_cell(current_adata)
	current_adata.raw = current_adata
	current_adata.layers["raw_count"] = current_adata.raw.X.copy()

	# Log transformation and scaling
	sc.pp.log1p(current_adata)
	sc.pp.scale(current_adata)

	# Load TF info which was made from mouse cell atlas dataset.
	base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()

	# Instantiate Oracle object
	oracle = co.Oracle()

	# In this notebook, we use the unscaled mRNA count for the input of Oracle object.
	current_adata.X = current_adata.layers["raw_count"].copy()

    # Instantiate Oracle object.
	oracle.import_anndata_as_raw_count(adata=current_adata,
	                               cluster_column_name="predicted.subclass",
	                               embedding_name="X_umap")

	# You can load TF info dataframe with the following code.
	oracle.import_TF_data(TF_info_matrix=base_GRN)

	# Perform PCA
	oracle.perform_PCA()

	# Select important PCs
	plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
	n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
	plt.axvline(n_comps, c="k")
	plt.show()
	print(n_comps)
	n_comps = min(n_comps, 50)

	n_cell = oracle.adata.shape[0]
	print(f"cell number is :{n_cell}")

	k = int(0.025*n_cell)
	print(f"Auto-selected k is :{k}")

	oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
	                  b_maxl=k*4, n_jobs=16)

	# Save oracle object.
	oracle.to_hdf5(f"/gpfs/home/asun/jinlab/deg_grn/ctxobj.{goi}_{sample}.{region}.{name}.subclass.celloracle.oracle")

	links = oracle.get_links(cluster_name_for_GRN_unit="predicted.subclass", alpha=10,
	                     verbose_level=0, n_jobs=16)

	links.to_hdf5(f"/gpfs/home/asun/jinlab/deg_grn/links.{goi}_{sample}.{region}.{name}.subclass.celloracle.links")

	links.filter_links()
	oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
	oracle.fit_GRN_for_simulation(alpha=10, 
	                          GRN_unit="cluster",
	                          use_cluster_specific_TFdict=True)

	oracle.to_hdf5(f"/gpfs/home/asun/jinlab/deg_grn/ctxobj.{goi}_{sample}.{region}.{name}.calculated.subclass.celloracle.oracle")