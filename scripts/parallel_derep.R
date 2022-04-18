# dereplicate and select O157 genomes in parallel

library(pdtools)
library(tidyverse)
library(tictoc)
library(furrr)
library(usethis)

if (future::supportsMulticore()){
  
  future::plan(multicore, workers=50)
  
} else {
  
  future::plan(multisession, workers=50)
  
}

options(future.globals.maxSize= 100048576000)



#######

pick_derep_sets <-
  function(meta, pan_PA, group=NULL, output_file, num_sets=25, desired_coverage=.99){
    
    if (!file.exists(output_file)){
      TIC <- tic()
      
      if (!is.null(group)){
        IDs <- meta %>% filter(GROUP == group) %>% pull(asm_acc)
        
        selection_PA <- pan_PA[,colnames(pan_PA) %in% IDs]
        selection_PA <- selection_PA[rowSums(selection_PA) > 0,]
        
        
      } else {
        
        selection_PA <- pan_PA
        selection_PA <- selection_PA[rowSums(selection_PA) > 0,]
        
      }
      
      derep_sets <-
        tibble(seed=seq(1:num_sets),
               # set90=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .90)),
               # set95=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .95)),
               selection_set=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = desired_coverage)),)
      
      saveRDS(derep_sets, output_file)
      print(group)
      TOC <- toc()
      print(TOC)
      
    } else {
      print('specified output file aready exists...returning it')
      derep_sets <- read_rds(output_file)
      
    }
    
    return(derep_sets)
    
    
  }



#######


use_directory('./output/')
# check and see if output sets are saved, if so, skip this
# if (!file.exists('./output/derep_sets_50_90.rds')){
#   
#   
# }

# read in presence absense (99% protein ID pangenome)
pan_PA <-
  read_tsv('./pan/WRITE/gene_presence_absence.Rtab') |> 
  column_to_rownames(var = 'Gene') |> 
  as.matrix()
meta <- read_tsv('./data/O157:H7_meta.tsv')

tic()

derep_sets <- 
  pick_derep_sets(meta = meta,
                  pan_PA = pan_PA,
                  num_sets = 100,
                  desired_coverage = .99, 
                  output_file = './output/derep_sets.rds')

TOC <- toc()

saveRDS(derep_sets, file = 'derep_sets.rds')

derep_sets <- read_rds('derep_sets.rds')

# gvt <- pan_mat_to_gene_vec_tibble(pan_PA)

# get_pangenome_representatives2(gvt)
# seed

# assembly accessions
# map(filtered_results_sets$set,1) 
# 
#
# set50_indexes <- map(derep_sets$set50,2)
# map(set50_indexes, length) %>% unlist() %>% hist()
# 
# set90_indexes <- map(derep_sets$set90,2)
# map(set90_indexes, length) %>% unlist() %>% hist()

# meta <- read_tsv('./data/O157:H7_meta.tsv')
# 
# meta <- 
#   meta %>%
#   extract_consensus_ag_species() %>%
#   right_join(meta) %>% 
#   extract_collection_agency() %>% 
#   right_join(meta)


# map(derep_sets$set90, 1)

# meta %>% filter(Year > 2018)

# # results_sets$set
# res_plot_dat <-
#   function(result, seed_num){
#     scores <- result[[2]]
#     # plot(1:length(scores), scores)
#     tibble(seed_num,num_genomes=1:length(scores), scores)
#   }
# 



# # 
# map2(.x = derep_sets$set50,.y = derep_sets$seed, res_plot_dat) %>%
#   bind_rows() %>%
#   ggplot(aes(x=num_genomes, y=scores, color=factor(seed_num))) + geom_point()



# 
# genome_set_frequency <-
#   map(derep_sets$set90, 1) %>%
#   purrr::reduce(c) %>%
#   tibble(genome_name=.) %>%
#   group_by(genome_name) %>%
#   summarise(asm_acc=unique(genome_name), 
#             num_derep_sets90=n(), .groups = 'drop') %>% 
#   select(-genome_name)
# 
# hist(genome_set_frequency$num_derep_sets90)
# 
# 
# genome_set_frequency %>% filter(num_derep_sets > 1)
# 
# 
# genome_set_frequency$asm_acc


# # 
# 
# filtered_results_sets2 %>%
#   mutate(genome_names=map(set75, 1),
#          scores=map(set75, 2)) %>%
#   select(-starts_with('set')) %>%
#   mutate(scores=map_chr(scores, ~paste(.x, collapse = ',')),
#          genome_names=map_chr(genome_names, ~paste(.x, collapse = ','))) %>%
#   write_tsv('./data/filtered_pan_derep75_50sets.tsv')
# 
# meta <- read_tsv('./data/O157:H7_meta.tsv')
# 
# meta_filt <- meta %>% right_join(genome_set_frequency)

# now want to check the coverage of this set



# 
# 
# set_PA <- pan_PA[rownames(pan_PA) %in% meta_filt$asm_acc,]
# dim(set_PA)
# 
# 
# set_PA_filt <- set_PA[,colSums(set_PA) > 0]
# dim(set_PA_filt)
# 
# ncol(set_PA_filt) / ncol(set_PA)

# 
# 
### need to extract CDC, FDA, FSIS(though i never saw these), etc from at least 2 columns



# 
# 
# meta_filt <- 
#   meta %>% 
#   left_join(genome_set_frequency) %>% 
#   # select(asm_acc, Year, bioproject_center,ag_match, collected_by,epi_type, STX, country, num_derep_sets90) %>% 
#   filter(!is.na(num_derep_sets90))

# consensus_sets %>% group_by()
# 
# meta_filt$collected_by
# meta_filt$bioproject_center
# 
# meta_filt %>% filter(grepl('FSIS', collected_by, ignore.case = TRUE))
# meta_filt$bioproject_center
# 
# 
# meta_filt %>% select(asm_acc,collection_agency,Year, PDS_acc, epi_type,ag_match)
# 
# 
# meta_filt$bioproject_center %>% unique()


######


#####

# 
# 
# 
# 
# 
# dim(pan_PA_filt)
# dim(pan_PA)
# 
# 
# library(parallelDist)
# jacc_dist <- parallelDist(pan_PA, method = 'binary')
# 
# 
# library(Rtsne)
# 
# 
# 
# 
# 
# 
# tsne_res <- Rtsne(is_distance=TRUE, X=jacc_dist)
# plot(tsne_res$Y)






