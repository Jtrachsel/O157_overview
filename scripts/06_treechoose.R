library(tidyverse)
library(pdtools)

source('./scripts/00_functions.R')
# metadata
meta <- read_tsv('data/O157:H7_meta.tsv')
# clusters
clusters <- read_tsv('./output/jaccard_clusters.tsv') %>% mutate(asm_acc=island_ID) %>% select(-island_ID)
# novelty
novelty <- read_tsv('./output/novelty_ranks.tsv')

dan_lin <- read_csv('output/DanLineages.csv', skip=1)
dan_lin %>% group_by(`LSPA type`) %>% tally() %>% arrange(desc(n))

dan_lin <- 
  dan_lin %>% 
  mutate(asm_acc=sub('data/genomes/O157/(.*).fna','\\1',Accession)) %>% 
  select(asm_acc, lineage)

meta <- meta %>% left_join(clusters) %>% left_join(novelty) %>% left_join(dan_lin)

tst <- 
  meta %>% 
  group_by(PDS_acc) %>% 
  summarise(prop_LI=sum(lineage %in% c('LI', 'LI/II'))/n(), 
            prop_LII=sum(lineage == 'LII')/n()) %>% 
  arrange(desc(prop_LII))

hist(tst$prop_LI)
hist(tst$prop_LII)

meta %>% write_tsv('output/05_meta.tsv')

meta %>% count(PDS_acc, lineage) %>%
  arrange(desc(n))


set.seed(7)
### this tree is too big
# tree_data <- 
#   meta %>%
#   filter(Year > 2017) %>%
#   filter(!is.na(PDS_acc)) %>% 
#   group_by(PDS_acc, Year, country, ag_match) %>% 
#   summarise(num_genomes=n(), 
#             representative=asm_acc[sample(1:n(), 1)]) %>%
#   arrange(desc(num_genomes))

###################
# reference genomes

# Nissle is too distantly related.  
# download.file('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/714/595/GCF_000714595.1_ASM71459v1/GCF_000714595.1_ASM71459v1_genomic.fna.gz', 
              # destfile='data/genomes/O157/Nissle.fna.gz')

# download.file('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/assembly_summary.txt',destfile = 'data/ecoli_refseq.txt')

# eco_ref <- read_tsv('./data/ecoli_refseq.txt', skip = 1)


# completes <- eco_ref %>% filter(assembly_level == 'Complete Genome')
# 
# library(pdtools)
# 
# potential_refs <- 
#   completes %>%
#   mutate(asm_acc=`# assembly_accession`) %>% 
#   make_download_urls('fna') %>%
#   make_dest_paths('fna','./data/references/') %>% 
#   download_genomes('fna')
# 
# potential_refs %>% filter(!fna_exists)


# fastANI -t 20 -q Sakai.fna --refList refseq_complete_ecoli.txt -o ecoli_vs_sakai.txt
# vs_sakai <- read_tsv('./data/references/ecoli_vs_sakai.txt', col_names = c('ref', 'query', 'ANI', 'v1', 'v2'))
# vs_sakai
# 
# vs_sakai <- vs_sakai %>%
#   mutate(query=sub('./(.*).fna','\\1', query)) 


# 
# 
# ref_types <- read_tsv('./data/references/ectyper/output.tsv')
# ref_types <- ref_types %>%
#   filter(Serotype != 'O157:H7') %>%
#   mutate(query=Name)
# 
# ref_types$query %in% vs_sakai$query

# vs_sakai <- vs_sakai %>% left_join(ref_types, by='query') %>% filter(!is.na(Serotype))

# vs_sakai %>% slice_head(n=2) %>% pull(query)

# 'GCF_013167975.1' = O55:H7 closest ANI relative to O157:H7s
# this one is maybe even a little too distant.
# system('cp ./data/references/GCF_013167975.1.fna ./data/genomes/O157/O55_H7.fna')
# hist(vs_sakai$ANI)

# when I plotted this tree I found that maybe 'GCA_011970025.1' would be a good outgroup
###


download.file('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz', 
              destfile = 'data/genomes/O157/Sakai.fna.gz')

download.file('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/168/095/GCF_013168095.1_ASM1316809v1/GCF_013168095.1_ASM1316809v1_genomic.fna.gz', 
              destfile='data/genomes/O157/86_24.fna.gz')

download.file('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/022/225/GCF_000022225.1_ASM2222v1/GCF_000022225.1_ASM2222v1_genomic.fna.gz',
              destfile = 'data/genomes/O157/TW14359.fna.gz')
# incomplete_genomes <- grep_vector(pattern_vector = tree_data$representative, list.files('./data/genomes/O157', full.names = T))
# complete_genomes <- grep('GCA_', list.files('./data/genomes/O157', full.names = T), invert = T, value = T)


# 
# build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes, 
#                              incomplete_genome_paths = incomplete_genomes) %>% 
#   write_tsv('overview_tree_ppanggolin.tsv', col_names = FALSE)
#########

# SET_SEED
set.seed(7)
tree_data_small <- 
  meta %>%
  filter(Year > 2017) %>%
  filter(!is.na(PDS_acc)) %>% 
  group_by(PDS_acc) %>% 
  summarise(num_genomes=n(), 
            representative=asm_acc[sample(1:n(), 1)]) %>%
  arrange(desc(num_genomes)) %>% 
  write_tsv('data/tree_reps.tsv')



incomplete_genomes <- grep_vector(pattern_vector = tree_data_small$representative, list.files('./data/genomes/O157', full.names = T))
complete_genomes <- grep('GCA_', list.files('./data/genomes/O157', full.names = T), invert = T, value = T)



build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes, 
                             incomplete_genome_paths = incomplete_genomes) %>% 
  write_tsv('small_overview_tree_ppanggolin.tsv', col_names = FALSE)



#########


