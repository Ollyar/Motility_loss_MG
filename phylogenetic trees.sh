#usr/bin/bash
#creates an alingment
muscle -align HFMG94VAA_RS03475_sequences.fasta -output HFMG94VAA_RS03475.aln

#create the tree
iqtree -s Virginia_consensus_strains_onlyLucy_1.fasta -B 1000 -T AUTO -m TEST --prefix iqtree_Virginia_consensus_strains_onlyLucy_1
