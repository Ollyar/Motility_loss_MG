#usr/bin/bash

#makes a blast database for all the protein sequences for the 78 MG strains
makeblastdb -dbtype prot -in Lucy_protein_db.faa -title MG_protine_sequences -out MG_protine_sequences -hash_index

#makes a blast database for the candidate genes we have extracted
makeblastdb -dbtype prot -in mutated_genes_with_protein_sequences.fasta -title Candidate_genes -out Candidate_genes -hash_index

blastp -query mutated_genes_with_protein_sequences.fasta -task blastp -db MG_protine_sequences -out Blastp_for_candidate_genes.tsv -evalue 1e-40 -outfmt 6 -show_gis -max_target_seqs 100

