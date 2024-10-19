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

def plot_scatter_with_labels(result, output_file, label):
    """
    Plots a scatter plot with gene labels for points where the absolute x-value is greater than 1.
    
    Parameters:
    result (pd.DataFrame): The dataframe containing 'Control_Log2_Fold_Change', 'logFC', and 'Gene' columns.
    output_file (str): The file path to save the output plot. Default is 'scatter_plot.jpg'.
    
    Returns:
    None: Displays and saves the scatter plot.
    """
    
    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=result, x='Temp_Log2_Fold_Change', y='logFC')
    
    # Add gene labels for significant points
    for i in range(len(result)):
        x_val = result['Temp_Log2_Fold_Change'][i]
        y_val = result['logFC'][i]
        if abs(x_val) > 1:
            plt.text(
                x_val, 
                y_val, 
                result['Gene'][i],
                fontsize=9,
                ha='right'  # horizontal alignment
            )
    
    # Set x-axis limits based on the maximum absolute value of 'Control_Log2_Fold_Change'
    max_x = max(abs(result['Temp_Log2_Fold_Change'].min()), abs(result['Temp_Log2_Fold_Change'].max()))
    plt.xlim(-max_x-1, max_x+1) 
    
    # Add labels, title, and reference lines
    plt.xlabel(f'{label} In Silico Log2 Fold Change')
    plt.ylabel('Functional Log2 Fold Change')
    plt.title(f'{label} vs Functional Log2 Fold Changes by Gene')
    
    plt.axvline(x=0, color='red', linestyle='--', linewidth=1)
    plt.axhline(y=0, color='red', linestyle='--', linewidth=1)
    
    # Save and display the plot
    plt.savefig(output_file, bbox_inches='tight', dpi=300)

goi_list= ["Tcf4"]  # Genes of interest
sample_list= ["1", "2", "3"]        # Sample numbers
region_list= ["CR", "Excit_L2_IT_ENTl", "Excit_L5_PT_CTX", "Excit_L5IT", "Excit_L5NP_CTX", "Excit_L6CT_CTX", "Excit_L6IT", "Excit_Upper", "Inhib_Sst"]
corresponding_region = ["CR", "L2 IT ENTl", "L5 PT CTX", "L5 IT CTX", "L5 NP CTX", "L6 CT CTX", "L6 IT CTX", "L2/3 IT CTX-1", "Sst"]
songpeng_regions = ["HPF_CR_Glut", "L2-3_IT_ENT_Glut", "", "L5_IT_CTX_Glut", "L5_NP_CTX_Glut", "L6_CT_CTX_Glut", "L6_IT_CTX_Glut", "", "", "Sst_Gaba"]
name_list = ["peturbseq", "ctrl", "Tbr1", "Foxg1", "Nr2f1", "Tcf4"]

# Loop through each list
for goi in goi_list:
    for sample in sample_list:
        for i in range(0, len(region_list)):
            region = region_list[i]

            for j in range(0, len(name_list)):
                name = name_list[j]
                azimuth_region = corresponding_region[i]
                
                file_path = f"/gpfs/home/asun/jinlab/deg_grn/calculated/ctxobj.{goi}_{sample}.{region}.{name}.calculated.subclass.celloracle.oracle"  # Replace with your file path
                if os.path.exists(file_path):
                    oracle = co.load_hdf5(f"/gpfs/home/asun/jinlab/deg_grn/calculated/ctxobj.{goi}_{sample}.{region}.{name}.calculated.subclass.celloracle.oracle")
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
                oracle_fc_df = create_fold_change_df(cluster_1, fold_change, log_fold_change, 'Temp')
                result = pd.merge(oracle_fc_df, DEG, on="Gene", how="inner")
                result = result[result['Temp_Log2_Fold_Change'] != 0]
                result.reset_index(drop=True, inplace=True)

                temp_fc = result['Temp_Log2_Fold_Change']
                functional_fc = result['logFC']  # LogFC from functional perturbations
                
                temp_class = [1 if fc > 0 else 0 for fc in temp_fc]
                functional_class = [1 if fc > 0 else 0 for fc in functional_fc]
                
                temp_tn, temp_fp, temp_fn, temp_tp = confusion_matrix(functional_class, temp_class).ravel()
                temp_fpr, temp_tpr, temp_thresholds = roc_curve(temp_class, functional_class)
                temp_roc_auc = auc(temp_fpr, temp_tpr)
                print(f"For {goi}_{sample}.{region}.{name}, AUC value of {temp_roc_auc:.2f}")
                
                if j == 0:
                    perturb_fpr = temp_fpr
                    perturb_tpr = temp_tpr
                    perturb_roc_auc = temp_roc_auc
                    plot_scatter_with_labels(result, f"/gpfs/home/asun/jinlab/graphs/scatterplot/{goi}_{sample}_{region}_{name}_scatter.jpg", "PerturbSeq")
                elif j == 1:
                    ctrl_fpr = temp_fpr
                    ctrl_tpr = temp_tpr
                    ctrl_roc_auc = temp_roc_auc
                    plot_scatter_with_labels(result, f"/gpfs/home/asun/jinlab/graphs/scatterplot/{goi}_{sample}_{region}_{name}_scatter.jpg", "Ctrl")
                else:
                    goi_fpr = temp_fpr
                    goi_tpr = temp_tpr
                    goi_roc_auc = temp_roc_auc
                    plot_scatter_with_labels(result, f"/gpfs/home/asun/jinlab/graphs/scatterplot/{goi}_{sample}_{region}_{name}_scatter.jpg", goi)

            plt.figure()
            songpeng_region = songpeng_regions[i]
            if songpeng_region != "":
                songpeng_oracle = co.load_hdf5(f"/gpfs/home/asun/jinlab/songpeng_grn/songpeng.{songpeng_region}.calculated.subclass.celloracle.oracle")
                
                # Enter perturbation conditions to simulate signal propagation after the perturbation.
                songpeng_oracle.simulate_shift(perturb_condition={goi: 0},
                                      GRN_unit="cluster",
                                      n_propagation=3)
                
                oracle_input = songpeng_oracle.adata.copy()
                oracle_count = songpeng_oracle.adata.copy()
                oracle_input.X = songpeng_oracle.adata.layers["simulation_input"].copy()
                oracle_count.X = songpeng_oracle.adata.layers["simulated_count"].copy()
                
                # Assuming 'adata' is your AnnData object and clusters are identified
                # 'cluster_key' is the key in adata.obs that contains cluster labels
                # 'cluster_1' and 'cluster_2' are the cluster labels you want to compare
                
                cluster_1_cells = oracle_count
                cluster_2_cells = oracle_input
                
                # Calculate mean expression for each gene in both clusters
                mean_expression_cluster_1 = cluster_1_cells.X.mean(axis=0)
                mean_expression_cluster_2 = cluster_2_cells.X.mean(axis=0)
                
                pseudocount = 1e-6 
                fold_change = (mean_expression_cluster_1 + pseudocount) / (mean_expression_cluster_2 + pseudocount)
                log_fold_change = np.log2((mean_expression_cluster_1 + pseudocount) / (mean_expression_cluster_2 + pseudocount))
                
                # Convert to a DataFrame for easy viewing
                genes = cluster_1_cells.var_names
                fold_change_df1 = pd.DataFrame({
                    'Gene': genes,
                    'Treatment_Fold_Change': fold_change,
                    'Treatment_Log2_Fold_Change': log_fold_change
                })
                
                # Sort by fold change or log fold change
                fold_change_df1 = fold_change_df1.sort_values(by='Treatment_Fold_Change', ascending=False)
                
                df = pd.read_csv(f"/gpfs/home/asun/jinlab/tsv/edgeR_LRT_with_sva.{goi}_{sample}.{region}.tsv", sep='\t')
                df = df.sort_values(by='logFC', ascending=False)
                df_sorted = df.iloc[df['logFC'].abs().argsort()]
                df_sorted.rename(columns={"Unnamed: 0": "Gene"}, inplace=True)
                merged_df = fold_change_df1.merge(df_sorted, how='left', left_on='Gene', right_on='Gene')
                merged_df.dropna()
                songpeng_fc = merged_df['Treatment_Log2_Fold_Change']
                functional_fc = merged_df['logFC']  # LogFC from functional perturbations
                
                songpeng_class = [1 if fc > 0 else 0 for fc in songpeng_fc]
                functional_class = [1 if fc > 0 else 0 for fc in functional_fc]

                songpeng_tn, songpeng_fp, songpeng_fn, songpeng_tp = confusion_matrix(functional_class, songpeng_class).ravel()
                
                songpeng_fpr, songpeng_tpr, songpeng_thresholds = roc_curve(songpeng_class, functional_class)
                songpeng_roc_auc = auc(songpeng_fpr, songpeng_tpr)
                print(f"For {goi}_{sample}.{region}, SONGPENG AUC value of {songpeng_roc_auc:.2f}")
                plt.plot(songpeng_fpr, songpeng_tpr, color='darkred', lw=2, label=f'Songpeng ROC curve (AUC = {songpeng_roc_auc:.2f})')

            plt.plot(ctrl_fpr, ctrl_tpr, color='darkorange', lw=2, label=f'Control ROC curve (AUC = {ctrl_roc_auc:.2f})')
            plt.plot(perturb_fpr, perturb_tpr, color='darkgreen', lw=2, label=f'PerturbSeq ROC curve (AUC = {perturb_roc_auc:.2f})')
            plt.plot(goi_fpr, goi_tpr, color='darkblue', lw=2, label=f'GOI ROC curve (AUC = {goi_roc_auc:.2f})')
            plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')  # Reference line for random guessing
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'AUROC for {goi} Perturbation in {region}')
            plt.legend(loc="lower right")
            plt.savefig(f"/gpfs/home/asun/jinlab/graphs/auroc/{goi}_{sample}_{region}_auroc.jpg", bbox_inches='tight', dpi=300)
                