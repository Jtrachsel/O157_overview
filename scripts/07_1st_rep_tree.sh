#!/bin/bash
set -e


ppanggolin annotate --fasta small_overview_tree_ppanggolin.tsv --cpu 40 -o smallpan -f
ppanggolin cluster -p smallpan/pangenome.h5 --cpu 40
ppanggolin graph -p smallpan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p smallpan/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p smallpan/pangenome.h5 -o smallpan/MSA --cpu 40 -f
ppanggolin rgp -p smallpan/pangenome.h5 -c 16
ppanggolin module -p smallpan/pangenome.h5 -c 16
ppanggolin spot -p smallpan/pangenome.h5 -c 16
ppanggolin write -o smallpan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p smallpan/pangenome.h5 -f
ppanggolin fasta -p ./smallpan/pangenome.h5 --output smallpan/REP_PROTS --prot_families all -f


raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n small_core -s ./smallpan/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7 -o GCA_011970025.1

