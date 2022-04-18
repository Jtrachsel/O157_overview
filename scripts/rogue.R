# tree_filt


library(ape)
library(Rogue)

rax_straps <- read.tree('RAxML_bootstrap.small_core')

# this was estimating 1 month for step 3....
rogue_output <- RogueTaxa(rax_straps,
                          bestTree = './RAxML_bipartitionsBranchLabels.small_core',
                          dropsetSize = 10)
rogue_output$taxon
