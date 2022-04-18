# cluster detection



suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(igraph))
library(digest)
library(parallelDist)
library(purrr)

# 
# overlap_dist <-  parallelDist::parallelDist(dat_mat, method='simpson') %>% 
#   write_rds('./gifrop_out/islands_overlap_dist.rds')

sim_mat <- 
  (1 - overlap_dist) %>%
  as.matrix() #%>%
# Matrix::Matrix(sparse = T)
rm(overlap_dist)

g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)
g <- delete_edges(g, E(g)[weight == 0])

print('writing overlap coef graph')
igraph::write.graph(g, file = './gifrop_out/overlap_coef_graph.dot', format = 'dot')
# find clusters in this network
print('Primary clustering, any islands sharing any number of genes will be in the same primary cluster')
print('two islands not sharing any gene content can be in the same primary cluster if they are linked 
          by and island that')
clust1 <- clusters(g)


print('secondary clustering')
print('pruning graph, removing edges with overlap coef of less than .5')
print('an overlap coefficient of .5 means that at least 1/2 of the genes in the smaller island are also present in the other island')
# print('any two islands with an overlap coef of at least .5 are in the same secondary cluster')

g <- delete_edges(g, E(g)[weight<scut])
clust2 <- cluster_louvain(g)

print('pruning graph, removing edges with overlap coef of less than XXX')


g <- delete_edges(g, E(g)[weight<tcut])
clust3 <- cluster_louvain(g)

print('constructing new graph with jaccard similarities')
print('This will allow the identification of extremely similar genomic islands')

# write dist objs for visualizations
jaccard_dist <- parallelDist::parallelDist(dat_mat, method = 'binary') %>% 
  write_rds('islands_jaccard_dist.rds')

sim_mat <- 
  (1 - jaccard_dist) %>% 
  as.matrix() #%>%
# Matrix::Matrix(sparse = T)
rm(jaccard_dist)

g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)
igraph::write.graph(g, file = './gifrop_out/jaccard_coef_graph.dot', format = 'dot')
print('removing edges representing jaccard similarities of less than .75')


g <- delete_edges(g, E(g)[weight<qcut])
clust4 <- cluster_louvain(g)

clust_info <- tibble(island_ID = names(membership(clust1)),
                     primary_cluster = membership(clust1),
                     secondary_cluster = membership(clust2),
                     tertiary_cluster = membership(clust3),
                     quat_cluster = membership(clust4))




#####

pan_PA <-
  read_tsv('./pan/gene_presence_absence.Rtab') |> 
  column_to_rownames(var = 'Gene') |> 
  as.matrix() |>
  t()

jaccard_dist <- parallelDist::parallelDist(pan_PA, method = 'binary')

write_rds(jaccard_dist, 'pan99_jaccard_dist.rds')


sim_mat <- 
  (1 - jaccard_dist) %>% 
  as.matrix()


g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)
# igraph::write.graph(g, file = './gifrop_out/jaccard_coef_graph.dot', format = 'dot')
# print('removing edges representing jaccard similarities of less than .75')

hist(sim_mat)

g <- delete_edges(g, E(g)[weight<.95])

louvain_clusts <- cluster_louvain(g)

clust_info <- 
  tibble(asm_acc=names(membership(louvain_clusts)), 
         cluster=membership(louvain_clusts))

factor(clust_info$cluster)


library(Rtsne)


TSNE <- Rtsne(jaccard_dist, is.dist=T)


TSNE
attributes(TSNE$Y)
attributes(jaccard_dist)$Labels

tsne_tib <- 
  tibble(asm_acc=attributes(jaccard_dist)$Labels, 
       TSNE1=TSNE$Y[,1],
       TSNE2=TSNE$Y[,2]
)


meta <- 
  read_tsv('data/O157:H7_meta.tsv') %>%
  left_join(clust_info) %>%
  left_join(tsne_tib) %>% 
  left_join(umap_tib)
library(tidyverse)
meta %>% ggplot(aes(x=TSNE1, y=TSNE2, color=fct_lump(factor(cluster), n = 8))) + geom_point()
meta %>% ggplot(aes(x=UMAP1, y=UMAP2, color=fct_lump(factor(cluster), n = 8))) + geom_point()

UMAP <- uwot::umap(jaccard_dist)

umap_tib <- 
  tibble(asm_acc=rownames(UMAP), 
       UMAP1=UMAP[,1], 
       UMAP2=UMAP[,2])




