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
	oracle = co.load_hdf5(f"/gpfs/home/asun/jinlab/songpeng_grn/songpeng.{group}.subclass.celloracle.oracle")
	links = co.load_hdf5(f"/gpfs/home/asun/jinlab/songpeng_grn/songpeng.{group}.subclass.celloracle.links")
	links.filter_links()
	oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
	oracle.fit_GRN_for_simulation(alpha=10, 
	                              GRN_unit="cluster",
	                              use_cluster_specific_TFdict=True)
	oracle.to_hdf5(f"/gpfs/home/asun/jinlab/songpeng_grn/songpeng.{group}.calculated.subclass.celloracle.oracle")