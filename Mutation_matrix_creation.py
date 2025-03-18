#%%
import vcfpy
import os
import pandas as pd
from BCBio import GFF

# Path to the directory containing the VCF files
vcf_dir = '/home/ollyar/data/sequences/Variant_results/SnpEff_results_virginia'
genomic_file = '/home/ollyar/data/sequences/genomic.gff'
metadata_file = '/home/ollyar/Metadata_genomes_lucy.csv'

# Initialize dictionaries to store gene positions and mutation data
gene_positions = {}
mutation_data = {}


# Step 1: Use BCBio.GFF to parse the genomic file and store gene and pseudogene positions
for record in GFF.parse(genomic_file):
        for feature in record.features:
            if feature.type in ['gene', 'pseudogene']:
                gene_name = feature.qualifiers.get('ID', [None])[0]
                if gene_name:
                    gene_name = gene_name.replace('gene-', '').replace('pseudogene-', '')
                    start_pos = int(feature.location.start + 1)  # Convert to 1-based position
                    end_pos = int(feature.location.end)
                    gene_positions[gene_name] = (start_pos, end_pos)

# Step 2: Loop through each VCF file and extract mutation data
for vcf_file in os.listdir(vcf_dir):
    if vcf_file.endswith('.vcf') or vcf_file.endswith('.vcf.gz'):
        full_vcf_path = os.path.join(vcf_dir, vcf_file)
        strain_name = os.path.basename(vcf_file).split('.')[0]
        mutation_data[strain_name] = []

        reader = vcfpy.Reader.from_path(full_vcf_path)

        for record in reader:
            mutation_position = record.POS
            effect_info = record.INFO.get('EFF', 'No Effect Info')

            # Filter for non-synonymous and stop codon mutations
            if effect_info != 'No Effect Info':
                for effect in effect_info:
                    if "NON_SYNONYMOUS" in effect or "STOP_GAINED" in effect or "STOP_LOST" in effect:
                        mutation_data[strain_name].append({
                            'Position': mutation_position,
                            'Effect': effect
                        })

# Step 3: Create a matrix for genes and strains
results_matrix = {gene: {strain: 0 for strain in mutation_data} for gene in gene_positions}

# Populate the results matrix
for strain, mutations in mutation_data.items():
    for mutation in mutations:
        mutation_position = mutation['Position']
        mutation_found = False

        for gene, (start_pos, end_pos) in gene_positions.items():
            if start_pos <= mutation_position <= end_pos:
                results_matrix[gene][strain] = 1  # Mark mutation presence in the gene
                mutation_found = True
                break  # Exit loop once match is found

        # Track mutations with no matching gene position
        if not mutation_found:
            print(f"No matching gene found for mutation at position {mutation_position} in strain {strain}")

# Step 4: Convert the results matrix to a DataFrame with genes on the x-axis
df = pd.DataFrame(results_matrix).T  # Transpose: genes as rows initially
df.fillna(0, inplace=True)
df = df.astype(int)

# Transpose the DataFrame to make genes the x-axis (columns) and strains the rows
df_transposed = df.T

# Save the transposed DataFrame to a CSV file
output_csv_transposed = '/home/ollyar/mutation_matrix.csv'
df_transposed.to_csv(output_csv_transposed)

# Output results
print(f"Data saved with genes on the x-axis to {output_csv_transposed}")


#%%
# Step 5: Load metadata and integrate sampling year
df_metadata = pd.read_csv(metadata_file)
df_metadata['strain'] = df_metadata['Folder.name'].str.extract(r'Sample_(\w+)$')
df_metadata.loc[df_metadata['strain'].isna(), 'strain'] = df_metadata['Folder.name']

strain_year_dict = dict(zip(df_metadata['strain'], df_metadata['year.sampling']))

# Add year column to the mutation matrix
df_transposed.insert(0, 'Year', df_transposed.index.map(strain_year_dict))

# Save the final DataFrame to a CSV file
output_csv_transposed = '/home/ollyar/mutation_matrix_with_years.csv'
df_transposed.to_csv(output_csv_transposed)

# Output results
print(f"Data saved with year column to {output_csv_transposed}")


