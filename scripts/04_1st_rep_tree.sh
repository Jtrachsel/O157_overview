#!/bin/bash
# must be run with "bash -i" because of 'conda activate'
set -e



# this is dan's python2 script.  grumble.
conda activate python2.7

LSPA6Long.sh output/tree_reps_lin.csv ./data/dan_lin/*fna

conda deactivate

#######

conda activate ppanggolin


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


# installed root_digger to find root of tree

#rd --msa ./smallpan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.small_core --prefix ROOT_DIG1 --threads 40
rd --msa ./smallpan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.small_core --prefix ROOT_DIG1 --threads 40 --exhaustive --early-stop


### TRY IQ TREE, has built in rooting options
# infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
#iqtree -p ALN_DIR --prefix concat -B 1000 -T AUTO

# infer the locus trees
#iqtree -S ALN_DIR --prefix loci -T AUTO

# compute concordance factors
#iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10