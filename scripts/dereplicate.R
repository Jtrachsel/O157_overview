library(pdtools)
library(tidyverse)

# install.packages("furrr")

# pangenome stuff?

# plan(multisession, workers = 20)

# options(future.globals.maxSize= 1048576000)
# 1000*1024^2 = 1048576000


pan_PA <-
  read_tsv('./pan/gene_presence_absence.Rtab') |> 
  column_to_rownames(var = 'Gene') |> 
  as.matrix() |>
  t()


# pan_PA[1:10, 1:10]

# pan_PA <- pan_PA[,colSums(pan_PA) > 2]


# reps <- get_gene_content_reps(pan_mat = pan_PA)
# 
# reps500 <- get_gene_content_reps(pan_mat = pan_PA,
#                               starting_set_size = 5000,
#                               desired_coverage = .98,
#                               set_size_step = 100,
#                               num_iters_per_size=500)
# 
# reps$value
# reps$set_size

# # pan_PA_filt <- pan_PA[,colSums(pan_PA) > 1]
# 
# singleton_genes <- pan_PA[,colSums(pan_PA) == 1]
# # singleton_genes <- singleton_genes[rowSums(singleton_genes) > 2, ]
# singleton genes from genomes with more than 1 singleton gene

# 
# # number of singleton genes per genome
# hist(rowSums(singleton_genes), breaks = 100)
# 
# # these genomes have no singleton genes in them
# no_singletons <- singleton_genes[rowSums(singleton_genes) == 0,]
# all(no_singletons == 0)
# dim(no_singletons)
# 
# # some singletons
# some_singletons <- singleton_genes[rowSums(singleton_genes) != 0,]
# rowSums(some_singletons) |> hist()

# 
# # rowSums(singleton_genes)
# # 
# meta <- read_tsv('./data/O157:H7_meta.tsv')
# 
# ag_species <- meta |> pdtools::extract_consensus_ag_species()
# 
# meta <- meta %>% left_join(ag_species)
# 
# meta <- rowSums(singleton_genes) |> 
#   enframe(name = 'asm_acc', value = 'num_singleton') |> 
#   right_join(meta)
# 
# 
# test <- 
#   singleton_genes |> 
#   as.data.frame() |> 
#   rownames_to_column(var='genome') |> 
#   gather(-genome, key = 'gene', value = 'presence') |> 
#   as_tibble()
# 
# 
# keep_these_genomes <- 
#   test |> 
#   group_by(genome) |> 
#   summarise(num_singletons=sum(presence)) |>
#   filter(num_singletons >= 5) |> pull(genome)
# 
# keep_these_singletons <- 
#   test |> filter(genome %in% keep_these_genomes) |> filter(presence == 1) |> 
#   pull(gene)
# 
# discard_these_singletons <- test |> 
#   filter(!(genome %in% keep_these_genomes)) |> filter(presence == 1) |> 
#   pull(gene)
# 

# rownames(pan_PA)

#keep these 
# 
# # get rid of these genes
# pan_PA_filt <- pan_PA[,(!(colnames(pan_PA) %in% discard_these_singletons))]
# 
# 
# 



# singleton_genes[rowSums(singleton_genes) > 4,] |>  colnames()
# # keep singletons that occur in genomes with more than 4 singletons
# # keep amr and virulence singletons
# # get rid of singletons that occur in genomes with less than 5 singletons
# 
# 
# meta |> filter(num_singleton > 1) |>  
#   # ggplot(aes(x=num_singleton, y=asm_stats_contig_n50)) + geom_hex()
#   # ggplot(aes(x=num_singleton, y=asm_stats_n_contig)) + geom_hex()
#   ggplot(aes(x=num_singleton, y=asm_stats_length_bp)) + geom_hex()
# 
# meta$asm_stats_contig_n50
# meta$asm_stats_length_bp
# meta$asm_stats_n_contig
# meta$assembly_method
# 
# meta |> ggplot(aes(x=))
# meta$Platform |> unique()
# 
# 
# num_sings_mod <- glm(family ='poisson' ,
#                      data = meta,
#                      formula=num_singleton~asm_stats_contig_n50+asm_stats_length_bp+asm_stats_n_contig+Platform)
# 
# 
# summary(num_sings_mod)

# number of genomes per singleton gene
# should all be == 1
# all(colSums(singleton_genes) ==1)
# TRUE


# pdtools::get_gene_content_reps()

# filtered_reps <- pdtools::get_pangenome_representatives(pan_PA_filt)

remove_strict_core <- function(pan_PA){
  
  pan_PA_strict_core_removed <- pan_PA[,colSums(pan_PA) != max(colSums(pan_PA))]
  return(pan_PA_strict_core_removed)
  
}

pan_PA_filt <- remove_strict_core(pan_PA)

collapse_exact_matches <- function(pan_PA){
  # collapse genomes with exactly the same gene content into one representative
  
  pan_PA_filt <- pan_PA[,colSums(pan_PA) != 1]
  
  hash_tib <- 
    tibble(hash=apply(pan_PA_filt, 1, digest::digest), 
           genome_name=names(hash), 
           num_genes=rowSums(pan_PA))
  hash_tib %>% 
    group_by(hash) %>% 
    arrange(desc(num_genes)) %>% 
    summarise(representative=genome_name[1],
              num_genomes=n(),
              all_genomes=paste(genome_name, collapse = ';')) %>% 
    arrange(desc(num_genomes))
  
}

# gene_vec_tibble <- pan_mat_to_gene_vec_tibble(pan_PA_filt)
# tictoc::tic()
# tst2 <- get_pangenome_representatives2(gene_vec_tibble, desired_coverage = .5)
# tictoc::toc()
# #34.133 sec elapsed
# 
# tictoc::tic()
# tst <- get_pangenome_representatives(pan_PA_filt, desired_coverage = .5)
# tictoc::toc()
# # 33.326 sec elapsed
# 
# obs <- object.size(pan_PA_filt)
# print(obs, units='auto')


filtered_results_sets <- 
  tibble(seed=seq(1:25), 
         set25=map(.x = seed, ~ get_pangenome_representatives(pan_mat = pan_PA_filt, SEED = .x, desired_coverage = .25)), 
         set50=map(.x = seed, ~ get_pangenome_representatives(pan_mat = pan_PA_filt, SEED = .x, desired_coverage = .50)), 
         set75=map(.x = seed, ~ get_pangenome_representatives(pan_mat = pan_PA_filt, SEED = .x, desired_coverage = .75)), 
         set85=map(.x = seed, ~ get_pangenome_representatives(pan_mat = pan_PA_filt, SEED = .x, desired_coverage = .85)), 
         set95=map(.x = seed, ~ get_pangenome_representatives(pan_mat = pan_PA_filt, SEED = .x, desired_coverage = .95)))


# gvt <- pan_mat_to_gene_vec_tibble(pan_PA)

# get_pangenome_representatives2(gvt)
# seed

# assembly accessions
# map(filtered_results_sets$set,1) 
# 
# # num genes
# set_indexes <- map(filtered_results_sets$set,2) 
# map(set_indexes, length) %>% unlist() %>% hist()
# 
# 
# # meta <- read_tsv('./data/North_America_stx_metadata.tsv')

# reps <- meta %>% filter(asm_acc %in% TST[[1]])

# 
# # results_sets$set
# res_plot_dat <- 
#   function(result, seed_num){
#   scores <- result[[2]]
#   # plot(1:length(scores), scores)
#   tibble(seed_num,num_genomes=1:length(scores), scores)
# }
# 
# # res_plot_dat(results_sets$set[[1]])
# 
# map2(.x = filtered_results_sets$set,.y = filtered_results_sets$seed, res_plot_dat) %>%
#   bind_rows() %>%
#   ggplot(aes(x=num_genomes, y=scores, color=factor(seed_num))) + geom_point()
# 
# genome_set_frequency <- 
#   map(filtered_results_sets$set, 1) %>% 
#  purrr::reduce(c) %>% 
#   tibble(genome_name=.) %>% 
#   group_by(genome_name) %>% 
#   tally()
# 
# genome_set_frequency %>% filter(n == 10)
# 
# genome_set_frequency$n %>% hist()
# 
# filtered_results_sets %>%
#   mutate(genome_names=map(set, 1), 
#          scores=map(set, 2)) %>%
#   select(-set) %>%
#   mutate(scores=map_chr(scores, ~paste(.x, collapse = ',')), 
#          genome_names=map_chr(genome_names, ~paste(.x, collapse = ','))) %>%
#   write_tsv('./data/filtered_pan_derep_95.tsv')
# 
# # meta <- read_tsv('./data/O157:H7_meta.tsv') 
# 
# meta_filt <- meta %>% filter(asm_acc %in% genome_set_frequency$genome_name)
# 
# 
# ### need to extract CDC, FDA, FSIS(though i never saw these), etc from at least 2 columns
# 
# meta_filt$collected_by
# meta_filt$bioproject_center
# 
# meta %>% filter(grepl('FSIS', collected_by, ignore.case = TRUE))
# meta$bioproject_center
# 
# 
# meta %>% select(asm_acc, PDS_acc, collected_by, bioproject_center, epi_type,ag_match)
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
