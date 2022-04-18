#######################
### this PDS summary section should really take into account the phylogenetic lineage
# move to after new tree calculation without 6 weird SNP clusters


library(tidyverse)
library(ggtree)
library(treeio)
library(pdtools)
library(ape)
source('scripts/00_functions.R')

meta <- read_tsv('./output/05_meta.tsv')



### read raxml tree ###
tr <- read.raxml('./RAxML_bipartitionsBranchLabels.small_pan2')

ggtree(tr, layout = 'fan') +# geom_nodelab(aes(label=node))+
geom_hilight(node=525) + geom_highlight(node=864, fill='green') + geom_highlight(node=686, fill='red')




FSIS_USDA <- 
  meta %>%
  filter(grepl('FSIS|USDA', collection_agency))

FSIS_PDSs <- FSIS_USDA %>% pull(PDS_acc) %>% unique()


PDS_summary <- 
  meta %>% 
  filter(!is.na(PDS_acc)) %>% 
  filter(Year > 2017) %>% 
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


tree_data <- 
  tibble(asm_acc=tr@phylo$tip.label) %>%
  left_join(meta) %>% 
  left_join(PDS_summary) %>% 
  mutate(LABEL=asm_acc)

tree_data$lineage

ggtr <- ggtree(tr) %<+% tree_data
# tree_data$LABEL

P_O157_1 <- 
  ggtr +
  geom_hilight(node=525, alpha=.25) +
  geom_highlight(node=864, fill='green', alpha=.25) +
  geom_highlight(node=686, fill='red', alpha=.25)+
  geom_tippoint(aes(fill=lineage, alpha=Year, size=total_isolates), shape=21) + 
  # geom_point2(aes(subset= FSIS == 'TRUE'),position = position_nudge(x = .00009, y = 0), color='red', size=1) + 
  geom_text2(aes(subset= !grepl('GCA',LABEL), label=LABEL), nudge_x = .00005, size=2) + 
  geom_cladelabel(node=525, label='clade_one') + 
  geom_cladelabel(node=864, label='clade_two') + 
  geom_cladelabel(node=686, label='clade_three')
  
P_O157_1

C1 <- extract.clade(tr@phylo, node=525)
C2 <- extract.clade(tr@phylo, node=864)
C3 <- extract.clade(tr@phylo, node=686)



tree_data <- 
  tree_data %>%
  mutate(tree_clade=
           case_when(
             asm_acc %in% C1$tip.label ~ 'clade_one', 
             asm_acc %in% C2$tip.label ~ 'clade_two', 
             asm_acc %in% C3$tip.label ~ 'clade_three',
             TRUE ~ 'other'
             
           ))


PDS_summary <- 
  tree_data %>%
  select(PDS_acc, tree_clade) %>%
  right_join(PDS_summary) %>% 
  arrange(desc(total_isolates)) %>% 
  write_tsv('output/09_PDS_summary.tsv')
  # mutate(all_ag_matches=fct_inorder(all_ag_matches))




meta <-
  meta %>%
  left_join(PDS_summary) %>%
  write_tsv('output/09_meta.tsv')


