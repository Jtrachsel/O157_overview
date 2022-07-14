#!/bin/bash
# must be run with "bash -i" because of 'conda activate'
set -e



# this is dan's python2 script.  grumble.
conda activate python2.7

LSPA6Long.sh Nathan_requests/tree_reps_lin.csv ./Nathan_requests/dan_lin/*fna

LSPA6Long.sh Nathan_requests_all/all_lineages.csv ./Nathan_requests_all/dan_lin/*fna


conda deactivate

#######

conda activate ppanggolin


ppanggolin annotate --fasta maybe_requestable_reps2022-06-22_ppanggolin.tsv --cpu 40 -o Nathan_requests_reps -f
ppanggolin cluster -p Nathan_requests_reps/pangenome.h5 --cpu 40
ppanggolin graph -p Nathan_requests_reps/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p Nathan_requests_reps/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p Nathan_requests_reps/pangenome.h5 -o Nathan_requests_reps/MSA --cpu 40 -f
ppanggolin msa --source protein --partition core --phylo -p Nathan_requests_reps/pangenome.h5 -o Nathan_requests_reps/MSA --cpu 40 -f
ppanggolin rgp -p Nathan_requests_reps/pangenome.h5 -c 16
ppanggolin module -p Nathan_requests_reps/pangenome.h5 -c 16
ppanggolin spot -p Nathan_requests_reps/pangenome.h5 -c 16
ppanggolin write -o Nathan_requests_reps/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p Nathan_requests_reps/pangenome.h5 -f
ppanggolin fasta -p ./Nathan_requests_reps/pangenome.h5 --output Nathan_requests_reps/REP_PROTS --prot_families all -f


raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n Nathan_requests_reps -s ./Nathan_requests_reps/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7


# installed root_digger to find root of tree

#rd --msa ./smallpan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.small_core --prefix ROOT_DIG1 --threads 40
rd --msa ./Nathan_requests_reps/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.Nathan_requests_reps --prefix Nathan_request_reps --threads 40 --exhaustive --early-stop

#### try iqtree ####


# To apply a codon model one should use the option -st CODON to tell IQ-TREE that the alignment contains protein coding sequences
# CODON11	The Bacterial, Archaeal and Plant Plastid Code
# ModelFinderPlus: -m MFP is the default behavior. Thus, this run is equivalent to iqtree -s example.phy.
# MSA dir: Nathan_requests_reps/MSA/msa_core_dna/
# -B 1000 UFboot,
#  A good practice is to use -T AUTO to determine the best number of cores:
# --seqtype STRING     BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)
  # -s DIR               Directory of alignment files

# need to screen the alignments before iqtree
# have to remove stop codons at least
# maybe run trimal 1st? remove any column with gaps?

# trim all columns containing gaps from the alignment


# codon model is a mess...
# try regular DNA
Rscript screen_alignment.R Nathan_requests_reps/MSA/msa_core_dna

iqtree --prefix NATE_REPS_IQ -T AUTO -B 1000 -s Nathan_requests_reps/MSA/msa_core_dna_screened/
iqtree --prefix NATE_REPS_IQ_S -T AUTO -B 1000 -S Nathan_requests_reps/MSA/msa_core_dna_screened/


#


trimal_all () {
  # $1 = input directory
  # #2 = output directory
  # mkdir -p $2
  echo $1
  
  for x in $(find "$1" -name "*.aln" )
  do
    BASE=$(basename -s .aln $x)
    out_file=$(printf '%s/%s_trim.aln' $2 $BASE)
    #echo $out_file
    trimal -nogaps -in $x -out "$out_file"
  done
  
}

# mkdir Nathan_requests_reps/MSA/msa_core_dna_trim/

# trimal_all Nathan_requests_reps/MSA/msa_core_dna Nathan_requests_reps/MSA/msa_core_dna_trim 
  

# Rscript remove_stop_codons.R Nathan_requests_reps/MSA/msa_core_dna_trim 

# iqtree --prefix NATE_REPS_IQ --seqtype CODON11 -T AUTO -B 1000 -s Nathan_requests_reps/MSA/msa_core_dna_no_stops/

# need to replace stop codons with ---  

# stop_codon = "(TAA|TGA|TAG)$")


### that was representatives for a core tree

# now do the full pangenome



ppanggolin annotate --fasta maybe_requestable_all2022-06-22_ppanggolin.tsv --cpu 40 -o Nathan_requests_all -f
ppanggolin cluster -p Nathan_requests_all/pangenome.h5 --cpu 40 --coverage .95 --identity .99
ppanggolin graph -p Nathan_requests_all/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p Nathan_requests_all/pangenome.h5
# ppanggolin msa --source dna --partition core --phylo -p Nathan_requests_all/pangenome.h5 -o Nathan_requests_all/MSA --cpu 40 -f
ppanggolin rgp -p Nathan_requests_all/pangenome.h5 -c 16
ppanggolin module -p Nathan_requests_all/pangenome.h5 -c 16
ppanggolin spot -p Nathan_requests_all/pangenome.h5 -c 16
ppanggolin write -o Nathan_requests_all/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p Nathan_requests_all/pangenome.h5 -f
ppanggolin fasta -p ./Nathan_requests_all/pangenome.h5 --output Nathan_requests_all/REP_PROTS --prot_families all -f
ppanggolin fasta -p ./Nathan_requests_all/pangenome.h5 --output MY_GENES --genes all


# annotate the protein families in the genomes
amrfinder -p Nathan_requests_all/REP_PROTS/all_protein_families.faa --plus -o Nathan_requests_all/amrfinder_rep_prots.tsv
psortb -i Nathan_requests_all/REP_PROTS/all_protein_families.faa -n --outdir Nathan_requests_all -o terse
# PSORTb version 3.0




### TRY IQ TREE, has built in rooting options
# infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
#iqtree -p ALN_DIR --prefix concat -B 1000 -T AUTO

# infer the locus trees
#iqtree -S ALN_DIR --prefix loci -T AUTO

# compute concordance factors
#iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10