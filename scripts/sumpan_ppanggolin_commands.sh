#!/bin/bash
set -e

# ppanggolin annotate --fasta overview_tree_ppanggolin.tsv --cpu 40 -o sumpan -f
ppanggolin cluster -p sumpan/pangenome.h5 --cpu 40
ppanggolin graph -p sumpan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p sumpan/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p sumpan/pangenome.h5 -o sumpan/MSA --cpu 40
ppanggolin rgp -p sumpan/pangenome.h5 -c 16
ppanggolin module -p sumpan/pangenome.h5 -c 16
ppanggolin spot -p sumpan/pangenome.h5 -c 16
ppanggolin write -o sumpan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p sumpan/pangenome.h5
ppanggolin fasta -p ./sumpan/pangenome.h5 --output sumpan/REP_PROTS --prot_families all
