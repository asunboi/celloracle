import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import anndata as ad
import seaborn as sns
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.metrics import roc_curve, auc

# Function to simulate shifts
def simulate_shift(oracle, goi, shift_value, GRN_unit="cluster", n_propagation=3):
    oracle.simulate_shift(perturb_condition={goi: shift_value}, 
                          GRN_unit=GRN_unit, 
                          n_propagation=n_propagation)
    return oracle

# Function to prepare oracle inputs and counts
def prepare_oracle_data(oracle):
    oracle_input = oracle.adata.copy()
    oracle_count = oracle.adata.copy()
    oracle_input.X = oracle.adata.layers["simulation_input"].copy()
    oracle_count.X = oracle.adata.layers["simulated_count"].copy()
    return oracle_input, oracle_count

# Function to calculate fold change and log2 fold change
def calculate_fold_changes(cluster_1_cells, cluster_2_cells, pseudocount=1e-6):
    mean_expression_cluster_1 = cluster_1_cells.X.mean(axis=0)
    mean_expression_cluster_2 = cluster_2_cells.X.mean(axis=0)
    fold_change = (mean_expression_cluster_1 + pseudocount) / (mean_expression_cluster_2 + pseudocount)
    log_fold_change = np.log2(fold_change)
    return fold_change, log_fold_change

# Function to create DataFrame for fold changes
def create_fold_change_df(cluster_1_cells, fold_change, log_fold_change, label):
    genes = cluster_1_cells.var_names
    fc_df = pd.DataFrame({
        'Gene': genes,
        f'{label}_Fold_Change': fold_change,
        f'{label}_Log2_Fold_Change': log_fold_change
    })
    return fc_df.sort_values(by=f'{label}_Fold_Change', ascending=False)

goi_list= ["Tcf4"]  # Genes of interest
sample_list= ["1", "2", "3"]        # Sample numbers
region_list= ["CR", "Excit_L2_IT_ENTl", "Excit_L5_PT_CTX", "Excit_L5IT", "Excit_L5NP_CTX", "Excit_L6CT_CTX", "Excit_L6IT", "Excit_Upper", "Inhib_Sst"]
corresponding_region = ["CR", "L2 IT ENTl", "L5 PT CTX", "L5 IT CTX", "L5 NP CTX", "L6 CT CTX", "L6 IT CTX", "L2/3 IT CTX-1", "Sst"]
name_list = ["peturbseq", "ctrl", "Tbr1", "Foxg1", "Nr2f1", "Tcf4"]

# Loop through each list
for goi in goi_list:
    for sample in sample_list:
        for i in range(0, len(region_list)):
                for name in name_list:
                    region = region_list[i]
                    azimuth_region = corresponding_region[i]
                    
                    file_path = f"/gpfs/home/asun/jinlab/deg_grn/ctxobj.{goi}_{sample}.{region}.{name}.calculated.subclass.celloracle.oracle"  # Replace with your file path
                    if os.path.exists(file_path):
                        oracle = co.load_hdf5(f"/gpfs/home/asun/jinlab/deg_grn/ctxobj.{goi}_{sample}.{region}.{name}.calculated.subclass.celloracle.oracle")
                    else:
                        print(f"Skipping {goi}_{sample}.{region}.{name}")
                        continue
                        
                    # Get Sean's DEG list
                    df = pd.read_csv(f"/gpfs/home/asun/jinlab/tsv/edgeR_LRT_with_sva.{goi}_{sample}.{region}.tsv", sep='\t')
                    df_sorted = df.iloc[df['logFC'].abs().argsort()]
                    df_sorted.set_index("Unnamed: 0")
                    rows_to_add = df_sorted.loc[df_sorted["Unnamed: 0"] == goi]
                    DEG = df_sorted.tail(3000)
                    DEG = pd.concat([DEG, rows_to_add], ignore_index=True)
                    DEG.rename(columns={'Unnamed: 0': 'Gene'}, inplace=True)
                        
                    oracle = simulate_shift(oracle, goi, 0)
                    oracle_input, oracle_count = prepare_oracle_data(oracle)
                    cluster_1 = oracle_count[oracle_count.obs['predicted.subclass'] == azimuth_region]
                    cluster_2 = oracle_input[oracle_input.obs['predicted.subclass'] == azimuth_region]
                    fold_change, log_fold_change = calculate_fold_changes(cluster_1, cluster_2)
                    oracle_fc_df = create_fold_change_df(cluster_1, fold_change, log_fold_change, 'PerturbSeq')
                    result = pd.merge(oracle_fc_df, DEG, on="Gene", how="inner")
                    result = result[result['PerturbSeq_Log2_Fold_Change'] != 0]
                    
                    perturb_fc = result['PerturbSeq_Log2_Fold_Change']
                    functional_fc = result['logFC']  # LogFC from functional perturbations
                    
                    perturb_class = [1 if fc > 0 else 0 for fc in perturb_fc]
                    functional_class = [1 if fc > 0 else 0 for fc in functional_fc]
                    
                    perturb_tn, perturb_fp, perturb_fn, perturb_tp = confusion_matrix(functional_class, perturb_class).ravel()
                    perturb_fpr, perturb_tpr, perturb_thresholds = roc_curve(perturb_class, functional_class)
                    perturb_roc_auc = auc(perturb_fpr, perturb_tpr)
                    print(f"For {goi}_{sample}.{region}.{name}, AUC value of {perturb_roc_auc:.2f}")