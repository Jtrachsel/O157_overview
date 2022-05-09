#!/bin/bash
set -e

ppanggolin annotate --fasta ppangg_file.tsv --cpu 30 -o pan -f
ppanggolin cluster -p pan/pangenome.h5 --cpu 30
ppanggolin graph -p pan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p pan/pangenome.h5
#ppanggolin msa --source dna --partition core --phylo -p pan/pangenome.h5 -o pan/MSA --cpu 40
ppanggolin rgp -p pan/pangenome.h5 -c 16
ppanggolin module -p pan/pangenome.h5 -c 16
ppanggolin spot -p pan/pangenome.h5 -c 16
ppanggolin write -o pan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p pan/pangenome.h5
ppanggolin fasta -p ./pan/pangenome.h5 --output pan/REP_PROTS --prot_families all
