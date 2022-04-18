#!/bin/bash
set -e


ppanggolin annotate --fasta second_overview_tree_ppanggolin.tsv --cpu 40 -o smallpan2 -f
ppanggolin cluster -p smallpan2/pangenome.h5 --cpu 40
ppanggolin graph -p smallpan2/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p smallpan2/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p smallpan2/pangenome.h5 -o smallpan2/MSA --cpu 40 -f
ppanggolin rgp -p smallpan2/pangenome.h5 -c 16
ppanggolin module -p smallpan2/pangenome.h5 -c 16
ppanggolin spot -p smallpan2/pangenome.h5 -c 16
ppanggolin write -o smallpan2/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p smallpan2/pangenome.h5 -f
ppanggolin fasta -p ./smallpan2/pangenome.h5 --output smallpan2/REP_PROTS --prot_families all -f


raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n small_pan2 -s ./smallpan2/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7 -o GCA_011970025.1

