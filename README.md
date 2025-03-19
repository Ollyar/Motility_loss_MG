# MG motility loss scripts

This repository contains scripts used to analyze MG genome and create figures for our dissertation. Raw data is not provided.

- [Mutation_matrix_creation.py](Mutation_matrix_creation.py) - Generate a mutation matrix displaying the presence or absence of mutations in all 789 genes across 78 strains. The matrix should indicate whether each gene in each strain contains at least one non-synonymous mutation, stop codon mutation, or start codon mutation.
- [Filter_faster_file_creation.py](Filter_faster_file_creation.py) - Creates a fasta file containing all genes mutated in non-motile strains but not in motile strains, corresponding protein sequences and specify which non-motile strains have these mutations.
- [Blast_database.sh](Blast_database.sh) - Outputs BLASTP results for candidate genes. 
- [Extracting_strain_gene_sequences.py](Extracting_strain_gene_sequences.py) - Creates a series of FASTA files, for each candidate gene, which contains the protine sequences for the homologus strain genes found using BLASTP.
- [phylogenetic trees.sh](phylogenetic trees.sh) - Creates an alingment for all strain genes, then run the alignment through IQ-TREE's model finder, and output the results.
- [Subplot_dis_graph.py](Subplot_dis_graph.py) - Creates Figure 1.

----------------------------------------------------------------
Acknowledgement: Thanks to Lily Smith and Poppy Daly for help with the code.
