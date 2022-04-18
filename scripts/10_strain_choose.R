library(tidyverse)
library(treeio)
library(ggtree)
library(readxl)
###
# read in NARMS data because we might be able to get these isolates more easily


meta <- read_tsv('./output/09_meta.tsv')
PDS_summary <- read_tsv('output/09_PDS_summary.tsv')


NARMS <- read_xlsx('NARMS_O157.xlsx', col_types = c('text'))%>%
  mutate(biosample_acc=`NCBI Accession Number`)


FSIS_USDA <- 
  meta %>%
  filter(grepl('FSIS|USDA', collection_agency))

FSIS_PDSs <- FSIS_USDA %>% pull(PDS_acc) %>% unique()

tr <- read.raxml('./RAxML_bipartitionsBranchLabels.small_pan2')


###

# picked biggest 3 SNP clusters from each of the 3 tree clades
# 9 SNP clusters under consideration



meta <- meta %>%
  mutate(NARMS=ifelse(biosample_acc %in% NARMS$biosample_acc, T, F), 
         FSIS=grepl('FSIS', isolate_identifiers))



complete_genomes <- 
  meta %>%
  filter(asm_level == 'Complete Genome') %>% 
  select(AMR_genotypes)



selected_SNPs <- 
  PDS_summary %>% 
  group_by(tree_clade) %>% 
  arrange(desc(total_isolates)) %>% 
  slice_head(n = 3) %>% 
  filter(grepl('clade_', tree_clade))


selected_SNPs %>% select(PDS_acc, total_isolates, hosts_except_human, USDA)
### show them on tree ###

# 
# interestingPDS <- 
#   meta %>% 
#   filter(!is.na(PDS_acc)) %>% 
#   filter(Year > 2017) %>% 
#   filter(grepl('clade',tree_clade)) %>% 
#   group_by(tree_clade, PDS_acc, Year) %>%
#   tally() %>% 
#   arrange(desc(n)) %>% 
#   # ungroup() %>% 
#   # group_by(Year) %>% 
#   mutate(tot_genomes=sum(n), 
#          prop_tot=n/sum(n) * 100) %>% 
#   filter(Year <= 2021) %>% ungroup() %>% 
#   group_by(PDS_acc, tree_clade) %>% 
#   summarise(biggest_year=Year[which.max(prop_tot)], 
#             tot_size=sum(n)) %>% 
#   filter(biggest_year == 2021) %>% 
#   arrange(desc(tot_size)) %>% 
#   ungroup() %>% group_by(tree_clade) %>% 
#   slice_head(n=2)
# 

# these two snp clusters maybe on the rise?
# biggest year by proportion was 2021 and they are decently big...
# 1 PDS000035159.319 clade_two          2021      491
# 2 PDS000035339.134 clade_one          2021      156

# 
# meta %>% 
#   filter(PDS_acc %in% interestingPDS$PDS_acc) %>% 
#   filter(!is.na(PDS_acc)) %>% 
#   filter(Year > 2017) %>% 
#   filter(grepl('clade',tree_clade)) %>% 
#   group_by(tree_clade, PDS_acc, Year) %>%
#   tally() %>% 
#   arrange(desc(n)) %>% 
#   # ungroup() %>% 
#   # group_by(Year) %>% 
#   mutate(tot_genomes=sum(n), 
#          prop_tot=n/sum(n) * 100) %>% 
#   filter(Year <= 2021) %>%
#   arrange(desc(tot_genomes)) %>% 
#   ggplot(aes(x=Year, y=prop_tot, size=n)) + 
#   geom_point() + 
#   geom_line(aes(group=PDS_acc), size=1) + 
#   facet_wrap(~tree_clade, scales = 'free')



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
  mutate(USDA=ifelse(PDS_acc %in% FSIS_PDSs, TRUE, FALSE), 
         SELECTED=ifelse(PDS_acc %in% selected_SNPs$PDS_acc, 'selected', NA))

tree_data <- 
  tibble(asm_acc=tr@phylo$tip.label) %>%
  left_join(meta) %>% 
  left_join(PDS_summary) %>% 
  mutate(LABEL=asm_acc)

# tree_data$lineage

ggtr <- ggtree(tr) %<+% tree_data
# tree_data$LABEL


ggtr +
  geom_hilight(node=525, alpha=.25) +
  geom_highlight(node=864, fill='green', alpha=.25) +
  geom_highlight(node=686, fill='red', alpha=.25)+
  geom_point2(aes(subset=PDS_acc %in% selected_SNPs$PDS_acc),shape=21, size=4, fill='purple') 



look <- meta %>%
  count(bioproject_center) %>% arrange(desc(n))


# some USDA sequenced isolates
'United State Department of Agriculture'
'USDA'
'U.S. Department of Agriculture'

meta %>% 
  mutate(possibility=
           case_when(
             grepl('USDA', bioproject_center) ~ TRUE, 
             grepl('Agriculture', bioproject_center) ~ TRUE, 
             NARMS ~ TRUE, 
             
                    
                  ))

'Enteric Diseases Laboratory Branch, Centers for Disease Control'


selected_isolates <- 
  meta %>% 
  filter(FSIS | NARMS)
  filter(Year > 2017) %>% 
  filter(PDS_acc %in% selected_SNPs$PDS_acc) %>% 
  mutate(novelty_score=ifelse(is.na(novelty_score), 0,novelty_score), 
         biosample_link=paste0('https://www.ncbi.nlm.nih.gov/biosample/', biosample_acc)) %>% 
  select(PDS_acc, strain,Year,ag_match,geo_loc_name, tree_clade, lineage,
         collection_agency, collected_by, bioproject_center, sra_center,
         RANK, biosample_link, asm_level, asm_acc, isolate_identifiers, novelty_score, FSIS, NARMS)


selected_isolates %>% count(NARMS)

#%>%
# select(asm_acc, tree_clade, lineage, strain, RANK, ag_match, Year, collection_agency, biosample_link)


# for each of the three tree clades I selected the top 3 SNP clusters by total isolates
# for each of the 9 considered SNP clusters I selected the top 2 novel isolates


# top two novel isolates for each group

novel_isolates <- 
  selected_isolates %>% 
  filter(grepl('clade_', tree_clade)) %>% 
  group_by(tree_clade, PDS_acc) %>%
  arrange(desc(novelty_score)) %>% 
  slice_head(n=2)# %>% 
# select(asm_acc, tree_clade, lineage, strain, RANK, ag_match, Year, collection_agency, biosample_link)

novel_isolates


# 
non_novel_isolates <- 
  selected_isolates %>% 
  filter(grepl('clade_', tree_clade)) %>% 
  group_by(tree_clade, PDS_acc) %>%
  arrange((novelty_score)) %>% 
  slice_head(n=2) #%>% 
# select(asm_acc, tree_clade, lineage, strain, RANK, ag_match, Year, collection_agency, biosample_link)


non_novel_isolates




selected_isolates %>% 
  filter(grepl('FSIS', isolate_identifiers)) %>%
  filter(asm_level == "Complete Genome") %>% 
  select(asm_acc, tree_clade, lineage, strain, RANK, ag_match, Year, collection_agency, biosample_link)

#####



most_novel_USDA_isolates <- 
  selected_isolates %>% 
  filter(grepl('FSIS', strain)) %>% 
  filter(grepl('clade_', tree_clade)) %>% 
  group_by(tree_clade, PDS_acc) %>%
  arrange(desc(novelty_score)) %>% 
  slice_head(n=2) %>% 
  select(asm_acc, tree_clade, lineage, strain, RANK, ag_match, Year, collection_agency, biosample_link)

most_novel_USDA_isolates


#####


# lineage assignment summary #
P_O157_2 <- 
  tree_data %>%
  select(PDS_acc, tree_clade) %>% 
  right_join(meta) %>% 
  filter(!is.na(tree_clade)) %>% 
  group_by(lineage, tree_clade) %>% 
  tally() %>% 
  filter(tree_clade != 'other') %>% 
  mutate(tree_clade=factor(tree_clade, levels = c('clade_one', 'clade_two', 'clade_three'))) %>% 
  ggplot(aes(x=tree_clade, y=n, fill=lineage)) + 
  geom_col() + 
  ylab('number of genomes')

# requesting isolates?
# step 1: select SNP clusters of interest
# step 2: within SNP clusters of interest, select isolates


host_order <- 
  PDS_summary %>% 
  filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three')) %>% 
  group_by(all_ag_matches) %>% summarise(tot_iso=sum(total_isolates)) %>% 
  arrange(tot_iso) %>% 
  pull(all_ag_matches)

PDS_summary <- 
  PDS_summary %>% 
  mutate(all_hosts=factor(all_ag_matches, levels = host_order), 
         tree_clade=factor(tree_clade, levels=c('clade_one', 'clade_two', 'clade_three'))) %>% 
  mutate(isos_per_year=total_isolates/num_years)


P_O157_3 <- 
  PDS_summary %>% 
  filter(!is.na(tree_clade)) %>% 
  # filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three')) %>% 
  ggplot(aes(x=total_isolates, y=all_hosts, fill=most_recent_year)) +
  geom_col(color='black') + 
  facet_wrap(~tree_clade, ncol=1, scales='free_y') + 
  scale_fill_viridis_c() + 
  ggtitle('SNP clusters by phylogenetic position and host', 
          'each rectangle represents one SNP cluster')
P_O157_3
####### now with only SNP clusters with isolates available from USDA


P_O157_4 <- 
  PDS_summary %>% 
  filter(USDA == 'TRUE') %>% 
  filter(!is.na(tree_clade)) %>% 
  ggplot(aes(x=total_isolates, y=all_hosts, fill=most_recent_year)) +
  geom_col(color='black') + 
  facet_wrap(~tree_clade, ncol=1, scales='free_y') + 
  scale_fill_viridis_c() + 
  ggtitle('SNP clusters by phylogenetic position and host, USDA isolates available', 
          'each rectangle represents one SNP cluster')
P_O157_4


O157_plots <- list(P_O157_1, P_O157_2, P_O157_3, P_O157_4)

O157_plots %>% write_rds('output/O157_plots.rds')

"/home/Julian.Trachsel/Documents/O157_overview/output/O157_plots.rds"

PDS_summary %>% 
  filter(!is.na(tree_clade)) %>% 
  ggplot(aes(x=log(total_isolates), log(isos_per_year), fill=all_hosts)) +
  geom_point(aes(size=most_recent_year),shape=21, color='black')


PDS_summary %>% arrange(desc(total_isolates))
# https://github.com/niemasd/TreeCluster

phylo_dist <- ape::cophenetic.phylo(tr@phylo) %>%
  as.data.frame() %>% 
  rownames_to_column(var='asm_acc') %>% 
  pivot_longer(cols=-asm_acc, names_to = 'to', values_to='phydist')

tst <- ape::cophenetic.phylo(tr@phylo) %>% as.dist()
HCLUST <- hclust(tst)
plot(HCLUST)
as.dist(tst)
hist(phylo_dist$phydist)

# geom_point2(aes(subset= Year == 2021),position = position_nudge(x = .00019, y = 0), color='blue', size=1)



PDS_summary %>% filter(USDA == 'TRUE')

PDS_summary %>% pull(hosts_except_human) %>% unique()

PDS_summary %>% ggplot(aes(x=total_isolates)) + geom_histogram()

#######



#######



PDS_summary %>%
  group_by(all_ag_matches) %>%
  summarise(tot_genomes=sum(total_isolates)) %>%
  arrange((tot_genomes)) %>% 
  mutate(all_ag_matches=fct_inorder(all_ag_matches)) %>% 
  ggplot(aes(y=all_ag_matches, x=tot_genomes)) + 
  geom_col()


######
meta %>% group_by(quat_cluster) %>% tally() %>% arrange(desc(n)) %>% 
  ggplot(aes(x=n)) + geom_histogram()


# THE NUMBER OF SNP CLUSTERS EACH QUAT CLUST IS FOUND IN

SNP_clusts_per_quat <- 
  meta %>%
  group_by(PDS_acc, quat_cluster) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(quat_cluster) %>% 
  tally()

SNP_clusts_per_quat %>% ggplot(aes(x=n)) + geom_histogram()





SNP_clusts_per_tert <- 
  meta %>%
  group_by(PDS_acc, tertiary_cluster) %>%
  tally() %>% 
  ungroup() %>% 
  group_by(tertiary_cluster) %>% 
  tally()

SNP_clusts_per_tert %>% ggplot(aes(x=n)) + geom_histogram()





SNP_clusts_per_sec <- 
  meta %>%
  group_by(PDS_acc, secondary_cluster) %>%
  tally() %>% 
  ungroup() %>% 
  group_by(secondary_cluster) %>% 
  tally()

SNP_clusts_per_sec %>% ggplot(aes(x=n)) + geom_histogram()

######


meta %>% 
  filter(!is.na(PDS_acc)) %>% 
  select(PDS_acc, ends_with('cluster')) %>% 
  arrange(quat_cluster)

meta %>%
  select(ends_with('cluster')) %>%
  unique() %>% 
  group_by(primary_cluster, secondary_cluster, tertiary_cluster, quat_cluster) %>%
  tally()


meta %>%
  select(ends_with('cluster')) %>%
  unique() %>% 
  group_by(primary_cluster, secondary_cluster, tertiary_cluster, quat_cluster) %>%
  tally()
