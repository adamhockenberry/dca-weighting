#!/bin/bash

for afa_file in /Users/adamhockenberry/Projects/DCA_weighting/Data/psicov150_aln_pdb/aln_fasta_max1k_bootstrap/*.fasta; do
    tree_file="$(basename "$afa_file" .fasta).newick"
    echo $tree_file
    /Users/adamhockenberry/workspace/FastTree/FastTree -lg -gamma -nosupport -quiet -out $tree_file $afa_file
    wait
done

