# info about PDS accessions
library(tidyverse)
library(ggtree)
library(treeio)
library(pdtools)

source('scripts/00_functions.R')

meta <- read_tsv('./output/O157:H7_meta.tsv')
tree_dat <- read_tsv('output/tree_reps.tsv') %>% 
  mutate(asm_acc=representative)


### lineages need work ###
tree_lin <- read_csv('output/tree_reps_lin.csv', skip = 1) %>% 
  mutate(asm_acc=sub('./data/dan_lin/(.*).fna','\\1',Accession))

tree_lin %>% 
  group_by(lineage, `LSPA type`) %>%
  tally() %>% 
  arrange(desc(n))

###


### read raxml tree ###
tr <- read.raxml('./RAxML_bipartitionsBranchLabels.small_core')

ggtree(tr)
# hist(tr@phylo$edge.length, breaks = 100)

# which(tr@phylo$edge.length > quantile(tr@phylo$edge.length, probs = .99))

phylo_dist <- ape::cophenetic.phylo(tr@phylo) %>%
  as.data.frame() %>% 
  rownames_to_column(var='asm_acc') %>% 
  pivot_longer(cols=-asm_acc, names_to = 'to', values_to='phydist')


# attempt to predict lineage based on phylo dists
to_lins <- tree_lin %>%
  transmute(to=asm_acc, to_lineage=lineage, to_LSPA=`LSPA type`)


# in this block the "to" column is the reference, known lineages to be matched against the unknowns
closest_known_lineages <- 
  tree_lin %>%
  left_join(phylo_dist) %>%  
  select(asm_acc, to, `LSPA type`, lineage, phydist) %>% 
  left_join(to_lins) %>% 
  filter(to_lineage != 'manually assign') %>% 
  filter(lineage == 'manually assign') %>% 
  group_by(asm_acc, `LSPA type`) %>% 
  summarise(min_phy=phydist[which.min(phydist)], 
            closest_known_lin=to_lineage[which.min(phydist)], 
            closest_LSPA=to_LSPA[which.min(phydist)], 
            closest_known_rep=to[which.min(phydist)])

closest_known_lineages %>% arrange((min_phy)) 

unknown_LSPA_votes <- 
  closest_known_lineages %>% 
  group_by(`LSPA type`, closest_LSPA, closest_known_lin) %>%
  tally() %>%
  arrange(desc(n)) %>% 
  ungroup() %>% 
  transmute(
    unknown_LSPA = `LSPA type`, 
    closest_known_LSPA=closest_LSPA, 
    closest_known_lin, 
    number_occurrences=n
  )




# only one LSPA type has a disagreement about what lineage it should be
lineage_conflict <- 
  closest_known_lineages %>% 
  group_by(`LSPA type`, closest_known_lin) %>%
  tally() %>%
  arrange(`LSPA type`) %>% 
  count(`LSPA type`) %>%
  filter(n !=1) %>% 
  pull(`LSPA type`)


# the only LSPA type with a conflict has only one vote for lin 1 and 29 votes for lin 2
all_conflicted <- closest_known_lineages %>% filter(`LSPA type` %in% lineage_conflict) 

# 221213 should probably be lineage 2

unknown_lineage_assignments <- 
  unknown_LSPA_votes %>%
  group_by(unknown_LSPA) %>% 
  summarise(pred_lineage=closest_known_lin[which.max(number_occurrences)])

# unknown_lineage_assignments
unknown_lin_vec <- unknown_lineage_assignments$pred_lineage
names(unknown_lin_vec) <- unknown_lineage_assignments$unknown_LSPA

known_lins <- 
  tree_lin %>% 
  filter(lineage != 'manually assign') %>% 
  select(`LSPA type`, lineage) %>% 
  unique()

known_lin_vec <- known_lins$lineage
names(known_lin_vec) <- known_lins$`LSPA type`

all_lin_vec <- c(unknown_lin_vec, known_lin_vec)

# this now has a new column LINEAGE that has the predicted lineage,
# determined the appropriate lineage to assign by voting on the closest known lineage
# as determined by phylogenetic distance from the ML tree
tree_lin <- 
  tree_lin %>% 
  mutate(LINEAGE=all_lin_vec[as.character(`LSPA type`)])

tree_lin %>% 
  select(asm_acc,`foldD-sfmA`, Z5935, yhcG, rtcB, rbsB, `arp-iclR`, `LSPA type`, lineage, LINEAGE) %>% 
  filter(lineage == 'manually assign') %>% 
  write_tsv('output/manual_lins_for_dan.tsv')

#####

phyclip_clusts <- read_tsv('phyclip_clusts/cluster_optimal_parameter_cs6_fdr0.3_gam4.5_sol0_f0_zero-branch-length-collapsed_RAxML_bestTree.txt') %>% 
  transmute(asm_acc=TAXA, 
            phyclip_clust=CLUSTER, 
            CLUST=fct_lump_n(f = as.factor(CLUSTER), n = 9))

# phyclip_clusts %>% group_by(CLUSTER) %>% tally() %>% arrange(desc(n))
###
tree_dat <-
  tree_dat %>%
  left_join(tree_lin) %>%
  filter(!is.na(lineage)) %>% 
  left_join(phyclip_clusts)




####### probably dont need this block
# used to identify leaves with very large total distances to others
# outliers?
cumulative_dists <- 
  phylo_dist %>%
  group_by(asm_acc) %>% 
  summarise(total_distance=sum(phydist)) %>% 
  arrange(desc(total_distance)) 



cumulative_dists %>%
  ggplot(aes(x=total_distance)) +
  geom_histogram()

weirdos <-  cumulative_dists %>% slice_head(n=15) %>% left_join(meta)

hist(cumulative_dists$total_distance, breaks = 100)
outlier_tips <- cumulative_dists %>% slice_head(n=6) %>% pull(asm_acc)

q99 <- quantile(cumulative_dists$total_distance, .99)


cumulative_dists %>% filter(total_distance > q99)

library(ape)
#######################
### this PDS summary section should really take into account the phylogenetic lineage
# move to after new tree calculation without 6 weird SNP clusters


FSIS_USDA <- 
  meta %>%
  filter(grepl('FSIS|USDA', collection_agency))

FSIS_PDSs <- FSIS_USDA %>% pull(PDS_acc) %>% unique()


PDS_summary <- 
  meta %>% 
  filter(!is.na(PDS_acc)) %>% 
  # group_by(PDS_acc, Year, country, ag_match, continent) %>%
  # tally() %>%
  # arrange(desc(n)) %>%
  # ungroup() %>% 
  group_by(PDS_acc) %>% 
  # nest()
  summarise(num_years=length(unique(Year)),
            total_isolates=n(),
            most_recent_year=unique(Year[which.max(Year)]), 
            # most_recent_year_num=sum(n[which.max(Year)]), 
            # previous_year_num = sum(n[which(Year == Year[which.max(Year)] - 1)]),
            all_ag_matches=paste(unique(sort(ag_match)), collapse = '_'),
            continents=paste(unique(sort(continent)), collapse = '_'),
            .groups='drop') %>%
  mutate(hosts_except_human=sub('Human', '', all_ag_matches), 
         hosts_except_human=sub('__','_',hosts_except_human), 
         hosts_except_human=sub('_$','',hosts_except_human), 
         hosts_except_human=sub('^_','',hosts_except_human)) %>% 
  arrange(desc(num_years)) %>% 
  mutate(USDA=ifelse(PDS_acc %in% FSIS_PDSs, 'TRUE', 'FALSE'))


tr@phylo
is.binary(tr@phylo)

tibble(edge_length=tr@phylo$edge.length) %>%
  ggplot(aes(x=edge_length)) + 
  geom_histogram(bins=100) +
  xlim(0,.000002)


all(tr@phylo$edge.length > 1e-6)
min(tr@phylo$edge.length)


# should run phyclip with 1.000001e-06 as the value to collapse 
###############

grep('GCA_',tree_data$asm_acc, invert = T) %>% length()
grep('GCA_',tree_data$asm_acc, invert = F) %>% length()

tree_data <- 
  tibble(asm_acc=tr@phylo$tip.label) %>%
  left_join(tree_dat) %>% 
  left_join(PDS_summary) %>% 
  mutate(LABEL=asm_acc, 
         LABEL2=ifelse(grepl('GCA', LABEL), '', LABEL))

tree_data$LINEAGE
tree_data$LABEL2

ggtr <- ggtree(tr) %<+% tree_data
tree_data$LABEL

library(ggrepel)




cbp1 <- c(LI="#009E73", `LI/II`="#E69F00", LII="#0072B2")

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# 
# 
# cbp1 <- c(None="#999999", Both="#E69F00", STX1="#56B4E9", STX2="#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




ggtr +
  geom_tippoint(aes(fill=lineage, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007)
# ggsave('output/lineage_tree_not_fixed.jpeg')

ggtr +
  geom_tippoint(aes(fill=LINEAGE, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007) + 
  scale_fill_manual(values = cbp1)+
  guides(fill = guide_legend(override.aes = list(size = 5)))
  # geom_point2(aes(subset= FSIS == 'TRUE'),position = position_nudge(x = .00009, y = 0), color='red', size=1) + 
  # geom_text2(aes(subset= !grepl('GCA',LABEL), label=LABEL), nudge_x = .00005, size=2) 
  # geom_point2(aes(subset= Year == 2021),position = position_nudge(x = .00019, y = 0), color='blue', size=1)

ggsave('output/lineages_tree.jpeg', bg = 'white' )


ggtr +
  geom_tippoint(aes(fill=CLUST, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007)

ggsave('output/phyclip_tree.jpeg' )


tree_data %>% filter(is.na(CLUST))

tree_data %>%
  group_by(CLUST, continents) %>%
  tally() %>%
  ungroup() %>% 
  arrange(desc(n)) %>% 
  group_by(CLUST) %>% 
  slice_head(n=1)


#######
ggtr +
  geom_tippoint(aes(fill=hosts_except_human, alpha=most_recent_year, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007)
#####


cp <- collapse(ggtr, node=c(1081))
cp <- collapse(cp, node=c(873))
cp <- collapse(cp, node=c(1436))

cp + 
  geom_point2(aes(subset=(node %in% c(1081,873, 1436 ))), size=5, shape=23, fill="purple")+
  geom_tippoint(aes(fill=LINEAGE, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007)+
  scale_fill_manual(values = cbp1)+
  guides(fill = guide_legend(override.aes = list(size = 5)))





ggsave('output/collapse_lineages_tree.jpeg', bg='white' )


cp + 
  geom_point2(aes(subset=(node %in% c(1081,873, 1436 ))), size=5, shape=23, fill="steelblue")+
  geom_tippoint(aes(fill=CLUST, alpha=most_recent_year, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007)

######

# log linear?


seniors <- array(data = c(911, 44, 538, 456, 3, 2, 43, 279), 
                 dim = c(2,2,2), 
                 dimnames = list("cigarette" = c("yes","no"),
                                 "marijuana" = c("yes","no"),
                                 "alcohol" = c("yes","no")))

# ggtr + geom_nodelab(aes(label=node))
# zoomClade(ggtr, node = 528)

#872 for north
#716 for mid?
#549 south?

# ggtr <- ggtree(tr, layout = 'circular')

# ggtr + geom_text(aes(label=node), hjust=.9 , size=3, color='blue') + 
  # geom_tiplab(size=2)

# extract.clade(tr@phylo, node = 1051)
# 863 may be for lineage 1?
# zoomClade(ggtr, node = 863)

# 706 for lin I/II?

# zoomClade(ggtr, node = 1051, xexpand = 2)  +geom_tiplab()
# ID main clades?

###
# 
# weird_clade <- ape::extract.clade(tr@phylo, node = 706)
# 
# remove_these <- weird_clade$tip.label
# tr@phylo$tip.label[!(tr@phylo$tip.label %in% remove_these)]


#### make new tree reps without strange outliers ###

# 
# new_tree_data <- 
#   tree_data %>%
#   filter(!(asm_acc %in% remove_these))
# 
# 
# incomp <- new_tree_data %>% filter(grepl('GCA', asm_acc))
# refs <- new_tree_data %>% filter(!grepl('GCA', asm_acc))
# 
# 
# incomplete_genomes <- grep_vector(pattern_vector = incomp$asm_acc, list.files('./data/genomes/O157', full.names = T))
# complete_genomes <- grep('GCA_', list.files('./data/genomes/O157', full.names = T), invert = T, value = T)
# 
# 
# 
# build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes, 
#                              incomplete_genome_paths = incomplete_genomes) %>% 
#   write_tsv('second_overview_tree_ppanggolin.tsv', col_names = FALSE)
# 

#############
