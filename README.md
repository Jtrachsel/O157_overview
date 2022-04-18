# O157_overview

Methods so far:  

1) Download metadata of E.coli in NCBI pathogens database  
2) Download genomes with both StxA and StxB  
3) use ectyper to serotype these genomes  
4) select only O157 genomes  
5) pangenome of all O157s  
6) cluster genomes by cloud partition  
7) select representatives and calculate core genome based phylogeny  
  a. One genome per SNP cluster  
  b. Reference genomes of interest  
  c. ppanggolin pangenome, uses mafft for core genome alignment  
  c. RAxML maximum likelihood GTRGAMMA model of concatenated core genome alignment
8) 6 SNP clusters with very long branch lengths, removing these and re-calc tree
9) classify genomes into lineage 'I', 'I/II', and 'II' with Daniel's script
10) tree for CRIS meeting?