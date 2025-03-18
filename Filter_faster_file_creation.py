#%%
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Path definitions
genbank_file = '/home/ollyar/data/sequences/GenomicSequenses.gbff'

# Define motile and non-motile strains
motile_strains = ['A1', 'F1', 'F4', 'A10', 'E11', 'E12']
non_motile_strains = ['B2', 'A9', 'D8', 'C3', 'B8']

# Load the mutation matrix from the previous step
mutation_matrix = pd.read_csv('/home/ollyar/mutation_matrix.csv', index_col=0)

# Identify genes mutated in motile strains
genes_mutated_in_motile = mutation_matrix[motile_strains].sum(axis=1) > 0

# Identify genes mutated in non-motile strains
genes_mutated_in_nonmotile = mutation_matrix[non_motile_strains].sum(axis=1) > 0

# Keep only genes mutated in non-motile strains but NOT in motile strains
filtered_genes = genes_mutated_in_nonmotile & ~genes_mutated_in_motile
filtered_matrix_nonmotile = mutation_matrix.loc[filtered_genes, non_motile_strains]

# Generate gene lists for each non-motile strain
mutated_genes_by_strain = {}
for strain in non_motile_strains:
    mutated_genes_by_strain[strain] = filtered_matrix_nonmotile.index[filtered_matrix_nonmotile[strain] == 1].tolist()

# Identify genes with mutations in any non-motile strain and in which strains they occur
mutated_genes_info = []
for gene in filtered_matrix_nonmotile.index:
    strains_with_mutation = [strain for strain in non_motile_strains if filtered_matrix_nonmotile.at[gene, strain] == 1]
    if strains_with_mutation:
        mutated_genes_info.append((gene, ', '.join(strains_with_mutation)))

# Extract protein sequences for each gene
gene_protein_sequences = {}
for record in SeqIO.parse(genbank_file, "genbank"):
    for feature in record.features:
        if feature.type == 'CDS':
            gene_id = feature.qualifiers.get('locus_tag', [None])[0]
            if gene_id and gene_id in filtered_matrix_nonmotile.index:
                protein_seq = feature.qualifiers.get('translation', ["No protein sequence available"])[0]
                gene_protein_sequences[gene_id] = protein_seq

# Create DataFrame to compile results
mutated_genes_df = pd.DataFrame(mutated_genes_info, columns=['Gene', 'Strains'])

# Add protein sequences, handling missing entries
mutated_genes_df['Protein Sequence'] = mutated_genes_df['Gene'].map(gene_protein_sequences).fillna("No protein sequence available")

# Output CSV file
output_csv = '/home/ollyar/mutated_genes_with_protein_sequences.csv'
mutated_genes_df.to_csv(output_csv, index=False)
print(f"Data saved to {output_csv}")

# Output FASTA file
output_fasta = '/home/ollyar/mutated_genes_with_protein_sequences.fasta'
fasta_records = []

for _, row in mutated_genes_df.iterrows():
    if row['Protein Sequence'] != "No protein sequence available":
        record = SeqRecord(Seq(row['Protein Sequence']), id=row['Gene'], description=f"Mutated in strains: {row['Strains']}")
        fasta_records.append(record)

SeqIO.write(fasta_records, output_fasta, "fasta")
print(f"FASTA file saved to {output_fasta}")
