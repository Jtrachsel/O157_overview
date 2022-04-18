

calc_ordinations <- function(DIST){
  # DIST <- parallelDist(x=t(pan_PA), method=method)

  TSNE <- Rtsne(DIST, is_distance = T)
  UMAP <- umap(X=DIST)
  MDS <- cmdscale(DIST)
  NMDS <- vegan::metaMDS(DIST)

  MDS <- MDS %>% as.data.frame() %>% rownames_to_column(var = 'asm_acc') %>%
    transmute(asm_acc, MDS1=V1, MDS2=V2)

  NMDS <- NMDS$points %>% as.data.frame() %>% rownames_to_column(var = 'asm_acc') %>%
    transmute(asm_acc, NMDS1=MDS1, NMDS2=MDS2)

  TSNE <- tibble(asm_acc=attributes(DIST)$Labels,
                 TSNE1=TSNE$Y[,1],
                 TSNE2=TSNE$Y[,2])

  UMAP <- UMAP %>% as.data.frame() %>%  rownames_to_column(var='asm_acc') %>%
    transmute(asm_acc, UMAP1=V1, UMAP2=V2)

  ord_coords <- MDS %>% left_join(NMDS) %>% left_join(TSNE) %>% left_join(UMAP)
  return(ord_coords)
}



mark_outliers <- function(DIST, outlier_prob=.99){

  mainmat <- as.matrix(DIST)
  dists_to_mine <- rowSums(mainmat)
  outliers <- dists_to_mine > quantile(dists_to_mine, probs = outlier_prob)
  print(paste0('detected ',sum(outliers), ' outliers', ' at ', outlier_prob ,' prob'))
  tibble(asm_acc=names(outliers),
         is_outlier=outliers)

}

#' Cluster genomes at 4 levels from a pangenome gene PA matrix
#' Needs parallelDist, igraph,
#'  genomes as rows and genes as columns?
#'
#' @param dat_mat
#' @param prefix
#' @param scut
#' @param tcut
#' @param qcut
#'
#' @return
#' @export
#'
#' @examples
cluster_genomes <-
  function(dat_mat,
           pcut=0,
           scut=.5,
           tcut=.75,
           qcut=.75,
           DIST_METHOD='simpson',
           output_directory=NULL,
           write_dist=TRUE){
    # browser()
    dist_filename <- paste0(output_directory, DIST_METHOD, '_dist.rds')
    graph_file_name <- paste0(output_directory, DIST_METHOD, '_graph.rds')

    if (!file.exists(graph_file_name)){
      print('calculating distances')

      # write dist objs for visualization
      # %>%


      if(!file.exists(dist_filename)){
        DIST <-  parallelDist::parallelDist(dat_mat, method=DIST_METHOD)
      } else {
        DIST <- read_rds(dist_filename)
      }

      if(write_dist){

        write_rds(DIST, dist_filename)
      }
      # write_rds('./gifrop_out/islands_overlap_dist.rds')

      print('converting to similarities')

      sim_mat <-
        (1 - DIST) %>%
        as.matrix() #%>%
      # Matrix::Matrix(sparse = T)
      rm(DIST)

      print('building graph')
      # maybe errors when no edges have weight 0?
      g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)

      print('writing graph')
      write_rds(g, graph_file_name)
      # igraph::write.graph(g, file = graph_file_name, format = 'edgelist')

    }  else {
      print(paste0('graph already exists, reading graph from', graph_file_name))
      g <- read_rds(graph_file_name)
    }


    if ( any(E(g)[weight < pcut]) ){
      g <- delete_edges(g, E(g)[weight < pcut])
    }

    # any connection = same cluster
    clust1 <- cluster_louvain(g)


    if ( any(E(g)[weight<scut]) ){
      g <- delete_edges(g, E(g)[weight<scut])
    }

    clust2 <- cluster_louvain(g)

    # print('pruning graph, removing edges with overlap coef of less than XXX')

    if ( any(E(g)[weight<tcut]) ){
      g <- delete_edges(g, E(g)[weight<tcut])
    }


    clust3 <- cluster_louvain(g)

    if ( any(E(g)[weight<qcut]) ){
      g <- delete_edges(g, E(g)[weight<qcut])
    }

    clust4 <- cluster_louvain(g)

    clust_info <- tibble(asm_acc = names(membership(clust1)),
                         primary_cluster = membership(clust1),
                         secondary_cluster = membership(clust2),
                         tertiary_cluster = membership(clust3),
                         quat_cluster = membership(clust4))

    return(clust_info)
  }


selection_orders <-
  function(VECTOR){
    tibble(genome_name=VECTOR,
           selected_order=1:length(VECTOR))
  }


# wrapper function to pick de-replication sets for specified NARMS serotype designations
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



calculate_novelty <-
  function(selection_set_results){
    # browser()
    genome_summary <-
      selection_set_results %>%
      mutate(genome_vectors=map(selection_set, 1),
             selection_orders=map(genome_vectors, selection_orders)) %>%
      select(selection_orders) %>%
      unnest(selection_orders) %>%
      filter(selected_order !=1) %>%
      group_by(genome_name) %>%
      summarise(mean_rank=mean(selected_order),
                median_rank=median(selected_order),
                number_selections=n(),
                best_rank=min(selected_order),
                worst_rank=max(selected_order)) %>%
      arrange(median_rank) %>%
      mutate(novelty_score1=number_selections*((1/(mean_rank + median_rank))),
             novelty_score2=(number_selections + (number_selections/(mean_rank + median_rank))),
             novelty_score3=((number_selections^2 / (mean_rank + median_rank))),
             novelty_score4=((number_selections / (median_rank))),
             novelty_score5=((number_selections^2 / (median_rank)^2)),
             novelty_score6=number_selections/max(number_selections) / (median_rank/n())) %>%

      filter(!(number_selections == 1 & best_rank == 1)) %>%
      transmute(asm_acc=genome_name,
                median_rank,
                number_selections,
                best_rank,
                worst_rank,
                novelty_score=novelty_score4,
                log_novelty=log(novelty_score)) %>%
      arrange(desc(novelty_score)) %>%
      ungroup() %>%
      mutate(RANK=1:n())

    return(genome_summary)

  }




res_plot_dat <-
  function(result, seed_num){
    scores <- result[[3]]
    # plot(1:length(scores), scores)
    tibble(seed_num,num_genomes=1:length(scores), scores)
  }



grep_vector <- function(pattern_vector, search_vector, INVERT=F){
  
  matches <- unique(grep(paste(pattern_vector,collapse="|"), 
                         search_vector,invert = INVERT, value=TRUE))
  return(matches)
  
}