import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import anndata as ad
import seaborn as sns

groups = ["HPF_CR_Glut", "L2-3_IT_ENT_Glut", "L5_IT_CTX_Glut", "L5_NP_CTX_Glut", "L6_CT_CTX_Glut", "L6_IT_CTX_Glut", "Sst_Gaba"]

for group in groups:
	# Load Songpeng adata and TFdict
	adata = sc.read_h5ad(f"/gpfs/home/asun/jinlab/songpeng_data/sa2.allen.subclass.{group}.ann.hdf5")
	tfi_all = co.motif_analysis.load_TFinfo(f"/gpfs/home/asun/jinlab/songpeng_data/sa2subclass.{group}.pdc.celloracle.tfinfo")
	TFdict = tfi_all.to_dictionary()

	adata.layers["raw_count"] = adata.X.copy()
	adata.layers["logCPM"] = adata.X.copy()
	sc.pp.scale(adata)
	sc.tl.pca(adata, svd_solver = "arpack")
	sc.pp.neighbors(adata, n_neighbors = 4, n_pcs = 50)
	sc.tl.umap(adata)

	oracle = co.Oracle()
	adata.X = adata.layers["logCPM"].copy()
	oracle.import_anndata_as_normalized_count(adata = adata,
	                                          cluster_column_name = "subclass",
	                                          embedding_name = 'X_umap')
	oracle.import_TF_data(TFdict = TFdict)
	oracle.perform_PCA()
	n_comps = np.where(np.diff(np.diff(
	    np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0]
	# plt.axvline(n_comps, c="k")
	# plt.show()
	n_comps = min(n_comps, 50)
	n_cell = oracle.adata.shape[0]
	k = max(int(0.025 * n_cell), 4)

	# KNN Imputation
	oracle.knn_imputation(n_pca_dims = n_comps, k = k, balanced = True,
	                      b_sight = k * 8, b_maxl = k*4, n_jobs=4)

	# get GRN
	links = oracle.get_links(cluster_name_for_GRN_unit = "subclass",
	                         alpha = 10, verbose_level=10)

	links.filter_links()
	oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
	oracle.fit_GRN_for_simulation(alpha=10, 
	                              GRN_unit="cluster",
	                              use_cluster_specific_TFdict=True)
	
	oracle.to_hdf5(f"/gpfs/home/asun/jinlab/songpeng_grn/songpeng.{group}.subclass.celloracle.oracle")
	links.to_hdf5(f"/gpfs/home/asun/jinlab/songpeng_grn/songpeng.{group}.subclass.celloracle.links")
