#%%
import pandas as pd
from Bio import SeqIO
import os

# Step 1: Load the CSV file into a DataFrame
file_path = "/home/ollyar/Blastp_for_candidate_genes_filtered.csv"
df = pd.read_csv(file_path)

# Step 2: Define the input FASTA file path
input_fasta_path = '/home/ollyar/data/sequences/Lucy_protein_db.faa'

# Step 3: Create output directory if it doesn't exist
output_dir = '/home/ollyar/data/trees/'
os.makedirs(output_dir, exist_ok=True)

# Step 4: Get a list of unique candidate genes
unique_genes = df['candidate_genes'].unique()

# Step 5: Process each candidate gene
for candidate_gene in unique_genes:
    # Filter the DataFrame for the current candidate gene and extract its strains
    strains_list = df[df['candidate_genes'] == candidate_gene]['strains'].tolist()

    # Define the output FASTA file path for the current candidate gene
    output_fasta_path = os.path.join(output_dir, f'{candidate_gene}_sequences.fasta')

    # Extract sequences for the strains and save them to a new FASTA file
    with open(output_fasta_path, 'w') as output_fasta:
        # Iterate over each sequence record in the input FASTA file
        for record in SeqIO.parse(input_fasta_path, 'fasta'):
            # Check if the current record's ID exactly matches any strain ID
            if record.id in strains_list:
                # Write the matching sequence to the output file
                SeqIO.write(record, output_fasta, 'fasta')

    print(f"Extracted sequences for {candidate_gene} saved to {output_fasta_path}")

print("Process completed for all candidate genes.")
