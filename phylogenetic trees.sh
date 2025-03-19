#usr/bin/bash
#creates an alingment
muscle -align HFMG94VAA_RS03475_sequences.fasta -output HFMG94VAA_RS03475.aln

#create the tree
iqtree -s HFMG94VAA_RS03475.aln -B 1000 -T AUTO -m TEST --prefix iqtree_HFMG94VAA_RS03475
