# info about PDS accessions
library(tidyverse)
library(ggtree)
library(treeio)
library(pdtools)

source('scripts/00_functions.R')

meta <- read_tsv('./output/05_meta.tsv')



### read raxml tree ###
tr <- read.raxml('./RAxML_bipartitionsBranchLabels.small_core')


# hist(tr@phylo$edge.length, breaks = 100)

# which(tr@phylo$edge.length > quantile(tr@phylo$edge.length, probs = .99))

phylo_dist <- ape::cophenetic.phylo(tr@phylo) %>%
  as.data.frame() %>% 
  rownames_to_column(var='asm_acc') %>% 
  pivot_longer(cols=-asm_acc, names_to = 'to', values_to='phydist')

# hist(phylo_dist$phydist)

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
  group_by(PDS_acc, Year, country, ag_match) %>%
  tally() %>%
  arrange(desc(n)) %>%
  ungroup() %>% 
  group_by(PDS_acc) %>% 
  # nest()
  summarise(num_years=length(unique(Year)),
            total_isolates=sum(n),
            most_recent_year=unique(Year[which.max(Year)]), 
            most_recent_year_num=sum(n[which.max(Year)]), 
            previous_year_num = sum(n[which(Year == Year[which.max(Year)] - 1)]),
            all_ag_matches=paste(unique(sort(ag_match)), collapse = '_'),
            .groups='drop') %>%
  mutate(hosts_except_human=sub('Human', '', all_ag_matches), 
         hosts_except_human=sub('__','_',hosts_except_human), 
         hosts_except_human=sub('_$','',hosts_except_human), 
         hosts_except_human=sub('^_','',hosts_except_human)) %>% 
  arrange(desc(num_years)) %>% 
  mutate(USDA=ifelse(PDS_acc %in% FSIS_PDSs, 'TRUE', 'FALSE'))





###############

tree_data <- 
  tibble(asm_acc=tr@phylo$tip.label) %>%
  left_join(meta) %>% 
  left_join(PDS_summary) %>% 
  mutate(LABEL=asm_acc)

tree_data$lineage

ggtr <- ggtree(tr) %<+% tree_data
tree_data$LABEL

ggtr +
  geom_tippoint(aes(fill=lineage, alpha=Year, size=total_isolates), shape=21) + 
  # geom_point2(aes(subset= FSIS == 'TRUE'),position = position_nudge(x = .00009, y = 0), color='red', size=1) + 
  geom_text2(aes(subset= !grepl('GCA',LABEL), label=LABEL), nudge_x = .00005, size=2) 
  # geom_point2(aes(subset= Year == 2021),position = position_nudge(x = .00019, y = 0), color='blue', size=1)

ggtr + geom_nodelab(aes(label=node))
# zoomClade(ggtr, node = 528)

#872 for north
#716 for mid?
#549 south?

ggtr <- ggtree(tr, layout = 'circular')

ggtr + geom_text(aes(label=node), hjust=.9 , size=3, color='blue') + 
  geom_tiplab(size=2)

# extract.clade(tr@phylo, node = 1051)
# 863 may be for lineage 1?
# zoomClade(ggtr, node = 863)

# 706 for lin I/II?

# zoomClade(ggtr, node = 1051, xexpand = 2)  +geom_tiplab()
# ID main clades?

###
# 
weird_clade <- ape::extract.clade(tr@phylo, node = 706)
# 
remove_these <- weird_clade$tip.label
tr@phylo$tip.label[!(tr@phylo$tip.label %in% remove_these)]


#### make new tree reps without strange outliers ###


new_tree_data <- 
  tree_data %>%
  filter(!(asm_acc %in% remove_these))


incomp <- new_tree_data %>% filter(grepl('GCA', asm_acc))
refs <- new_tree_data %>% filter(!grepl('GCA', asm_acc))


incomplete_genomes <- grep_vector(pattern_vector = incomp$asm_acc, list.files('./data/genomes/O157', full.names = T))
complete_genomes <- grep('GCA_', list.files('./data/genomes/O157', full.names = T), invert = T, value = T)



build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes, 
                             incomplete_genome_paths = incomplete_genomes) %>% 
  write_tsv('second_overview_tree_ppanggolin.tsv', col_names = FALSE)


#############
