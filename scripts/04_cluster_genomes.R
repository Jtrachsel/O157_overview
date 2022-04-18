library(tidyverse)
library(igraph)
library(parallelDist)
###

###########
source('./scripts/00_functions.R')


###
meta <- read_tsv('./data/O157:H7_meta.tsv')

meta %>%
  group_by(PDS_acc, Year, ag_match, country) %>%
  tally() %>% arrange(desc(n))

pan_PA <-
  read_tsv('./pan/WRITE/gene_presence_absence.Rtab') |> 
  column_to_rownames(var = 'Gene') |> 
  as.matrix()

# this takes forever, should feed in column types?
gpa <- read_csv('./pan/WRITE/matrix.csv')

table(gpa$`Non-unique Gene name`)

clouds <- gpa %>% filter(`Non-unique Gene name` == 'cloud') %>% pull(Gene)

gpa_cloud <- pan_PA[rownames(pan_PA) %in% clouds,]

# not all genomes have cloud genes
gpa_cloud <- gpa_cloud[,colSums(gpa_cloud) > 0]

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
