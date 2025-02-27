import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import pandas as pd

goi = "Tbr1"

# Read data from the text file
with open(f'jinlab/{goi}_AUC.txt', 'r') as file:
    data = file.readlines()  # Read all lines into a list
# Initialize an empty list to store the parsed values
parsed_data = []

# Loop through each line in the data
for line in data:
    if line.startswith("For"):
        # Split the line into parts
        parts = line.split(", AUC value of ")
        
        # Extract the gene, region, subset information
        gene_region_subset = parts[0].split(" ")[1]  # Get the gene_region_subset part
        auc = float(parts[1])  # Extract the AUC value
        
        # Further split gene_region_subset by '_'
        gene, region_subset = gene_region_subset.split('_', 1)
        
        # Split region_subset by the '.' delimiter
        region, subset = region_subset.split('.', 1)

        brain, last = subset.split('.', 1)

        # Append the values to the list
        parsed_data.append([gene, region, brain, last, auc])

# Create a DataFrame from the parsed data
df = pd.DataFrame(parsed_data, columns=["Gene", "Replicate", "Region", "Subset", "AUC"])

# Pivot the DataFrame to get the maximum AUC for each gene, replicate, region, and subset
pivot_df = df.pivot_table(index=['Gene', 'Replicate', 'Region'], 
                          columns='Subset', 
                          values='AUC', 
                          aggfunc='max').reset_index()

# Determine colors based on conditions
def assign_color(row):
    if row['peturbseq'] > row['ctrl']:
        return 'blue'
    elif row['ctrl'] > row['peturbseq']:
        return 'red'
    else:
        return 'black'

# Apply color assignment
pivot_df['Color'] = pivot_df.apply(assign_color, axis=1)
pivot_df['max'] = pivot_df[['peturbseq', 'ctrl', goi]].max(axis=1)


# Prepare the heatmap color data
heatmap_data = pivot_df.pivot_table(index=['Replicate', 'Gene'], 
                                    columns='Region', 
                                    values='Color', 
                                    aggfunc='first')

# Prepare the heatmap AUC data for annotation
auc_data = pivot_df.pivot_table(index=['Replicate', 'Gene'], 
                                columns='Region', 
                                values='max',  # Use peturbseq as the representative AUC values
                                aggfunc='first')

# Plotting the heatmap
plt.figure(figsize=(10, 6))
ax = sns.heatmap(heatmap_data.applymap(lambda x: 2 if x == 'blue' else (0 if x == 'red' else 1)), 
                 cmap='RdBu', 
                 cbar=False, 
                 linewidths=.5, 
                 linecolor='gray', 
                 annot=auc_data,  # Annotate heatmap with AUC values
                 fmt='.2f',  # Format AUC values to two decimal points
                 annot_kws={'size': 10, 'weight': 'bold'})  # Customize annotation size

# Customize the plot
plt.title(f'{goi} AUC Heatmap')
plt.xlabel('Region')
plt.ylabel('Replicate and Gene')
plt.xticks(rotation=45)
plt.yticks(rotation=0)

# Create a legend
blue_patch = mpatches.Patch(color='navy', label='Perturbseq')
red_patch = mpatches.Patch(color='crimson', label='Ctrl')
black_patch = mpatches.Patch(color='white', label=goi)

# Add the legend to the plot
plt.legend(handles=[blue_patch, red_patch, black_patch], 
           title="Conditions", 
           loc='upper right', 
           bbox_to_anchor=(1.2, 1))  # Adjust position as needed


#plt.show()
plt.savefig(f"jinlab/graphs/{goi}_AUC_Heatmap.jpg")