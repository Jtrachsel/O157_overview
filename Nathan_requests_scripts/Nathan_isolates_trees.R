library(tidyverse)
library(ggtree)

maybe_requestable <- read_tsv('Nathan_requests/maybe_requestable_metadata.tsv')


tree_data <- read_tsv('Nathan_requests/maybe_requestable_tree_reps_2022-06-21.tsv') %>% 
  transmute(asm_acc=representative, PDS_acc)
# tree_data <- read_tsv('Nathan_requests/maybe_requestable_tree_reps_2022-06-27.tsv')
# tree_data <- read_tsv('Nathan_requests/maybe_requestable_tree_reps_2022-06-28.tsv')



tr <- read.tree('Nathan_request_reps.rooted.tree')
tr2 <- read.tree('NATE_REPS_IQ.contree')

ggtr <- ggtree(tr)
ggtr2 <- ggtree(tr2)
ggtr
# ggtr2



# all lineages
all_lins <- read_csv('Nathan_requests_all/all_lineages.csv', skip=1)

all_lins <- 
  all_lins %>%
  transmute(asm_acc=sub('./Nathan_requests_all/dan_lin/(.*).fna','\\1',Accession), 
            individual_lineage=lineage)

maybe_requestable <- 
  maybe_requestable %>%
  left_join(all_lins)

## representative tree lineages
lins <- read_csv('Nathan_requests/tree_reps_lin.csv', skip = 1)

lins <- lins %>%
  mutate(asm_acc=sub('./Nathan_requests/dan_lin/(.*).fna','\\1',Accession))


phylo_dist <- ape::cophenetic.phylo(tr) %>%
  as.data.frame() %>% 
  rownames_to_column(var='asm_acc') %>% 
  pivot_longer(cols=-asm_acc, names_to = 'to', values_to='phydist')


# attempt to predict lineage based on phylo dists
to_lins <- lins %>%
  transmute(to=asm_acc, to_lineage=lineage, to_LSPA=`LSPA type`)


# in this block the "to" column is the reference, known lineages to be matched against the unknowns
closest_known_lineages <- 
  lins %>%
  left_join(phylo_dist) %>%  
  select(asm_acc, to, `LSPA type`, lineage, phydist) %>% 
  left_join(to_lins) %>% 
  filter(to_lineage != 'manually assign') %>% 
  filter(lineage == 'manually assign') %>% 
  group_by(asm_acc, `LSPA type`) %>% 
  summarise(min_phy=phydist[which.min(phydist)], 
            closest_known_lin=to_lineage[which.min(phydist)], 
            closest_LSPA=to_LSPA[which.min(phydist)], 
            closest_known_rep=to[which.min(phydist)],
            .groups = 'drop')

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


# No LSPA type conflicts
lineage_conflict <- 
  closest_known_lineages %>% 
  group_by(`LSPA type`, closest_known_lin) %>%
  tally() %>%
  arrange(`LSPA type`) %>% 
  count(`LSPA type`) %>%
  filter(n !=1) %>% 
  pull(`LSPA type`)


# None
all_conflicted <- closest_known_lineages %>% filter(`LSPA type` %in% lineage_conflict) 

unknown_lineage_assignments <- 
  unknown_LSPA_votes %>%
  group_by(unknown_LSPA) %>% 
  summarise(pred_lineage=closest_known_lin[which.max(number_occurrences)])

# unknown_lineage_assignments
unknown_lin_vec <- unknown_lineage_assignments$pred_lineage
names(unknown_lin_vec) <- unknown_lineage_assignments$unknown_LSPA

known_lins <- 
  lins %>% 
  filter(lineage != 'manually assign') %>% 
  select(`LSPA type`, lineage) %>% 
  unique()

known_lin_vec <- known_lins$lineage
names(known_lin_vec) <- known_lins$`LSPA type`

all_lin_vec <- c(unknown_lin_vec, known_lin_vec)

# this now has a new column LINEAGE that has the predicted lineage,
# determined the appropriate lineage to assign by voting on the closest known lineage
# as determined by phylogenetic distance from the ML tree
lins <- 
  lins %>% 
  mutate(LINEAGE=all_lin_vec[as.character(`LSPA type`)]) %>% 
  select(asm_acc,`foldD-sfmA`, Z5935, yhcG, rtcB,
         rbsB, `arp-iclR`, `LSPA type`, lineage, LINEAGE)



lins %>% 
  filter(lineage == 'manually assign') %>% 
  write_tsv('Nathan_requests/manual_lineage_assignments.tsv')

tree_lins <- 
  lins %>%
  left_join(maybe_requestable) %>%
  select(asm_acc, PDS_acc, LINEAGE)

# these are the associations between lineage and SNP cluster from the 
# representatives tree

PDS_lins <-
  tree_lins %>% 
  transmute(PDS_acc, PDS_lineage=LINEAGE) %>% 
  filter(!is.na(PDS_lineage)) 


maybe_requestable <- maybe_requestable %>% left_join(PDS_lins)

maybe_requestable %>% write_tsv('Nathan_requests/maybe_requestable_metadata.tsv')

## ADD STATES ###

PDS_summary <- 
  maybe_requestable %>% 
  filter(!is.na(PDS_acc)) %>% 
  group_by(PDS_acc) %>% 
  summarise(num_years=length(unique(Year)),
            total_isolates=n(),
            most_recent_year=unique(Year[which.max(Year)]), 
            all_ag_matches=paste(unique(sort(ag_match)), collapse = '_'),
            States=paste(unique(sort(State)), collapse = '_'),
            countries=paste(unique(sort(country)), collapse = '_'),
            collection_agencies=paste(unique(sort(collection_agency)), collapse = '_'),
            lineages=paste(unique(sort(individual_lineage)), collapse = '_'),
            .groups='drop') %>%
  mutate(hosts_except_human=sub('Human', '', all_ag_matches), 
         hosts_except_human=sub('__','_',hosts_except_human), 
         hosts_except_human=sub('_$','',hosts_except_human), 
         hosts_except_human=sub('^_','',hosts_except_human)) %>% 
  arrange(desc(num_years)) 

write_tsv(PDS_summary, 'Nathan_requests/PDS_summary.tsv')

tree_data <- 
  tibble(asm_acc=tr$tip.label) %>%
  left_join(tree_data) %>%
  left_join(PDS_summary) %>% 
  left_join(tree_lins) %>% 
  mutate(LABEL=asm_acc, 
         LABEL2=ifelse(grepl('GCA', LABEL), '', LABEL))

tree_data$lineages
tree_data$LABEL2

ggtr <- ggtree(tr) %<+% tree_data

# to put lineage I isolates at the top of the tree
ggtr <- flip(ggtr, 369, 469) 

# ggtr + geom_tippoint(aes(fill=LINEAGE, size=total_isolates), shape=21)

library(ggrepel)
cbp1 <- c(LI="#009E73", `LI/II`="#E69F00", LII="#0072B2")

# lineages

ggtr +
  geom_tippoint(aes(fill=LINEAGE, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007) + 
  scale_fill_manual(values = cbp1)+
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave('Nathan_requests/lineage_tree.jpeg')
# hosts

# PDS accession tip labels

ggtr +
  geom_tiplab(aes(label=PDS_acc),size=2, shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00015) + 
  scale_fill_manual(values = cbp1)+
  guides(fill = guide_legend(override.aes = list(size = 5)))


ggsave('Nathan_requests/lineage_tree_PDS_labs.jpeg', height = 15, width = 7, units = 'in')


ggtr +
  geom_tippoint(aes(fill=all_ag_matches, size=total_isolates), shape=21) + 
  geom_text_repel(aes(label=LABEL2), max.overlaps = 10000, size=3, nudge_x = .00007) + 
  # scale_fill_manual(values = cbp1)+
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave('Nathan_requests/hosts_tree.jpeg')

# plot tree


