library(tidyverse)
library(igraph)
library(parallelDist)
###

###########
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
    
    clust_info <- tibble(island_ID = names(membership(clust1)),
                         primary_cluster = membership(clust1),
                         secondary_cluster = membership(clust2),
                         tertiary_cluster = membership(clust3),
                         quat_cluster = membership(clust4))
    
    return(clust_info)
  }



###
meta <- read_tsv('./data/O157:H7_meta.tsv')

meta %>%
  group_by(PDS_acc, Year, ag_match, country) %>%
  tally() %>% arrange(desc(n))

pan_PA <-
  read_tsv('./pan/WRITE/gene_presence_absence.Rtab') |> 
  column_to_rownames(var = 'Gene') |> 
  as.matrix()

gpa <- read_csv('./pan/WRITE/matrix.csv')

unique(gpa$`Non-unique Gene name`)
table(gpa$`Non-unique Gene name`)

# colnames(meta)

clouds <- gpa %>% filter(`Non-unique Gene name` == 'cloud') %>% pull(Gene)


gpa_cloud <- pan_PA[rownames(pan_PA) %in% clouds,]

# not all genomes have cloud genes
gpa_cloud <- gpa_cloud[,colSums(gpa_cloud) > 0]

# gpa_cloud %>% as.data.frame() %>% rownames

simpson_clusters <- cluster_genomes(t(gpa_cloud),
                             DIST_METHOD = 'simpson',
                             output_directory = './output/',
                             pcut = 0,
                             scut = .25,
                             tcut = .50, 
                             qcut=.75)

write_tsv(simpson_clusters, './output/simpson_clusters.tsv')


jaccard_clusters <- cluster_genomes(t(gpa_cloud),
                                    DIST_METHOD = 'binary',
                                    output_directory = './output/',
                                    pcut = 0,
                                    scut = .25,
                                    tcut = .50,
                                    qcut = .75)

write_tsv(jaccard_clusters, './output/jaccard_clusters.tsv')

jaccard_clusters %>% 
  group_by(primary_cluster, secondary_cluster) %>%
  tally() %>%
  arrange(desc(n)) %>% 
  filter(n >2)
