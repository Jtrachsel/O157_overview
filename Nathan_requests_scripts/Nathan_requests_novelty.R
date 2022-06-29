library(pdtools)
library(tidyverse)
library(tictoc)
library(furrr)
library(usethis)

pick_derep_sets_jaccard <-
  function(DIST, output_file, num_sets=25, MAX_GENOMES){
    
    if (!file.exists(output_file)){
      # TIC <- tic()
      
      # selection_PA <- pan_PA
      # selection_PA <- selection_PA[rowSums(selection_PA) > 0,]
      
      
      
      derep_sets <-
        tibble(seed=seq(1:num_sets),
               # set90=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .90)),
               # set95=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .95)),
               selection_set=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives_jaccard(pan_dist = DIST, SEED = .x, max_genomes = MAX_GENOMES)),)
      
      saveRDS(derep_sets, output_file)
      # print(group)
      # TOC <- toc()
      # print(TOC)
      
    } else {
      print('specified output file aready exists...returning it')
      derep_sets <- read_rds(output_file)
      
    }
    
    return(derep_sets)
    
    
  }


if (future::supportsMulticore()){
  
  future::plan(multicore, workers=25)
  
} else {
  
  future::plan(multisession, workers=25)
  
}

options(future.globals.maxSize= 100048576000)




amrfinder <- read_tsv('Nathan_requests_all/amrfinder_rep_prots.tsv')


amrfinder %>%
  group_by(`Gene symbol`) %>%
  tally() %>% 
  arrange(desc(n)) %>% 
  filter(grepl('stx', `Gene symbol`))



gpa <- read_tsv('Nathan_requests_all/WRITE/gene_presence_absence.Rtab')

pan_mat <- gpa %>% column_to_rownames(var = 'Gene') %>% as.matrix()


## these cover 99% of all the genes in the pangenome
## a 99% ID 95% coverage pangenome
## defaults are 80% ID 80% coverage (AA data)
greedy_sets <- 
  pdtools::pick_derep_sets(pan_mat, output_file = 'Nathan_requests_all/greedy_sets.rds')

greedy_novelty <- calculate_novelty(greedy_sets)


######
if(!file.exists('Nathan_requests_all/DIST.rds')){
  print('dist not found, generating....')
  pan_dist <- parallelDist::parallelDist(base::t(pan_mat), method = 'binary')
  write_rds(x = pan_dist, file = 'Nathan_requests_all/DIST.rds')
} else {
  pan_dist <- read_rds('Nathan_requests_all/DIST.rds')
}


jaccard_sets <- 
  pick_derep_sets_jaccard(DIST = pan_dist, MAX_GENOMES = 700, output_file = 'Nathan_requests_all/jaccard_sets.rds')

jaccard_novelty <- calculate_novelty(jaccard_sets)


## membrane novelty ##

psort_res <- 
  read_tsv('Nathan_requests_all/20220624140437_psortb_gramneg.txt')
  
psort_res$Localization %>% unique()

membrane_prots <- psort_res %>% 
  filter(!(Localization %in% c('Unknown', 'Cytoplasmic'))) %>% 
  pull(SeqID)

membrane_mat <- pan_mat[rownames(pan_mat) %in% membrane_prots,]

colSums(membrane_mat) %>% hist(breaks=1000)
rowSums(membrane_mat) %>% hist(breaks=1000)


## removes proteins found in all genomes

membrane_mat <- membrane_mat[rowSums(membrane_mat) != ncol(membrane_mat),]
ncol(membrane_mat)
rowSums(membrane_mat) %>% hist(breaks=100)



membrane_mat[1:5, 1:5]

membrane_greedy_sets <- 
  pdtools::pick_derep_sets(membrane_mat, output_file = 'Nathan_requests_all/membrane_greedy_sets.rds')

membrane_greedy_novelty <- calculate_novelty(membrane_greedy_sets)


#



###


ranked_genomes <- 
  greedy_novelty %>%
      transmute(asm_acc, greedy_RANK=RANK) %>% 
  full_join(
    membrane_greedy_novelty %>% transmute(asm_acc, membrane_greedy_RANK=RANK)
  ) %>% 
  arrange((greedy_RANK))




ranked_genomes %>% ggplot(aes(x=membrane_greedy_RANK, y=greedy_RANK))+geom_point()

####

meta <-
  read_tsv('./Nathan_requests/maybe_requestable_metadata.tsv') %>% 
  left_join(ranked_genomes) %>% 
  write_tsv('./Nathan_requests/maybe_requestable_metadata.tsv')

# meta %>% 
#   select(asm_acc, ag_match, Year, country, State, epi_type, isolation_source,
#          collection_agency, isolate_identifiers, NARMS, collected_by,
#          outbreak) %>% 
#   left_join(ranked_genomes) %>% 
#   ggplot(aes(x=ag_match, y=jaccard_RANK)) + geom_point() 
# 

meta %>% 
  select(asm_acc, ag_match, Year, country, State,
         collection_agency, isolate_identifiers, NARMS, collected_by,
         epi_type, outbreak) %>% 
  left_join(ranked_genomes) %>%
  filter(is.na(membrane_greedy_RANK) & is.na(greedy_RANK))


