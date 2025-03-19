#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the metadata file to get strain years
metadata_file = "/home/ollyar/Metadata_genomes_lucy.csv"
df_metadata = pd.read_csv(metadata_file)
df_metadata['strain'] = df_metadata['Folder.name'].str.extract(r'Sample_(\w+)$')
df_metadata.loc[df_metadata['strain'].isna(), 'strain'] = df_metadata['Folder.name']
strain_year_dict = dict(zip(df_metadata['strain'], df_metadata['year.sampling']))

# Load the mutation matrix
mutation_file = "/home/ollyar/mutation_matrix.csv"
df = pd.read_csv(mutation_file)

# Extract gene names and create numerical indices
genes = df.columns[1:]
gene_indices = np.arange(1, len(genes) + 1)  # Numbering genes from 1 to 789

# Pre/Post-2007 strain mutation counts
pre_2007_counts = {gene: 0 for gene in genes}
post_2007_counts = {gene: 0 for gene in genes}

total_pre_2007 = 0
total_post_2007 = 0

# Classify strains based on year and count mutations
for index, row in df.iterrows():
    strain_name = row.iloc[0]  # Get strain name from first column
    if strain_name in strain_year_dict:
        year = strain_year_dict[strain_name]
        if year < 2007:
            total_pre_2007 += 1
            for gene in genes:
                if row[gene] == 1:
                    pre_2007_counts[gene] += 1
        else:
            total_post_2007 += 1
            for gene in genes:
                if row[gene] == 1:
                    post_2007_counts[gene] += 1

# Convert mutation counts to proportions
pre_2007_proportions = [count / total_pre_2007 if total_pre_2007 > 0 else 0 for count in pre_2007_counts.values()]
post_2007_proportions = [count / total_post_2007 if total_post_2007 > 0 else 0 for count in post_2007_counts.values()]

# Convert dictionaries to lists for plotting
genes = list(pre_2007_counts.keys())  # Get gene names
gene_indices = np.arange(1, len(genes) + 1)  # Generate indices from 1 to 789

# Create figure with 2 subplots
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(20, 10), sharex=True)

# Plot Pre-2007 Proportion Bar Graph
axes[0].bar(gene_indices, pre_2007_proportions, color='blue', label='Pre-2007 Strains')
axes[0].set_ylabel("Proportion Strains with Mutations", fontsize=17)

# Plot Post-2007 Proportion Bar Graph
axes[1].bar(gene_indices, post_2007_proportions, color='red', label='Post-2007 Strains')
axes[1].set_xlabel("Gene Index (Ordered 1-789)", fontsize=18)
axes[1].set_ylabel("Proportion Strains with Mutations", fontsize=17)

# Improve readability: Show every 10th gene index for clarity
plt.xticks(ticks=gene_indices[::10], labels=gene_indices[::10], rotation=90, fontsize=10)

# Add a single legend below the graph
fig.legend(["Pre-2007 Strains", "Post-2007 Strains"], loc='lower center', ncol=2, fontsize=14)

# Adjust layout to fit everything properly
plt.tight_layout(rect=[0, 0.05, 1, 1])  # Give space for the legend

# Save the figure
plt.savefig("/home/ollyar/mutation_proportions_pre_post_2007_bars.png", dpi=300)

# Show the figure
plt.show()

# Extract gene names and create numerical indices
genes = df.columns[1:]  # Extract gene names from the columns
gene_indices = list(range(1, len(genes) + 1))  # Generate indices from 1 to 789

# Create a DataFrame mapping Gene Index to Gene Name
gene_table = pd.DataFrame({"Gene Index": gene_indices, "Gene Name": genes})

# Save table as a CSV file for inclusion in the appendix
table_path = "/home/ollyar/gene_index_table.csv"
gene_table.to_csv(table_path, index=False)

print(f"Table saved as: {table_path}")