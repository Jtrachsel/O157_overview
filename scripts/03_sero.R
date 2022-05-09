library(tidyverse)
library(cowplot)
library(Biostrings)
library(pdtools)

# ran ectyper on all downloaded genomes


ectyper_colnames <- 
  c('Name', 'Species', 'O-type',
  'H-type', 'Serotype', 'QC', 
  'Evidence', 'GeneScores', 'AlleleKeys', 'GeneIdentities(%)', 
  'GeneCoverages(%)', 'GeneContigNames', 'GeneRanges', 'GeneLengths', 
  'Database', 'Warnings')

# output of ectyper

sero <- 
  read_tsv('./data/all_sero.tsv', col_names = ectyper_colnames) %>% 
  mutate(asm_acc=sub('(GCA_[0-9]+\\.[0-9]).*','\\1',Name)) 
# because I downloaded and serotyped some genomes the old way i needed this mutate

# LOOK <- sero %>% filter(grepl('H7', `H-type`))

SERO <- sero %>% select(asm_acc, Serotype)

# output of download script
meta <-
  read_tsv('./output/all_stx_metadata.tsv') %>% 
  left_join(SERO)

meta <- meta %>% filter(!is.na(PDS_acc))

O157_PDS <- meta %>% filter(Serotype == 'O157:H7') %>% pull(PDS_acc) %>% unique()


meta %>%
  filter(PDS_acc %in% O157_PDS) %>% 
  count(PDS_acc, Serotype) %>%
  ungroup() %>% 
  arrange(desc(n)) %>% 
  count(PDS_acc) %>% arrange(desc(n))


# meta %>% filter(PDS_acc == 'PDS000004368.66') %>% select(Serotype)
# meta %>% filter(!(asm_acc %in% sero$Name)) %>% pull(asm_acc)
# sero %>% filter(!(Name %in% meta$asm_acc)) %>% pull(Name)
# sero %>% filter(!(asm_acc %in% meta$asm_acc)) %>% pull(asm_acc)

sero %>% 
  select(asm_acc, Serotype) %>% 
  left_join(meta) %>%
  mutate(serotype=fct_lump((Serotype), n = 10)) %>%
  group_by(serotype, STX) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>%
  mutate(serotype2=fct_inorder(serotype)) %>%
  ggplot(aes(y=serotype, x=n)) + 
  geom_col(aes(fill=STX)) + 
  ggtitle('serotypes of stx containing E.coli') + 
  theme_cowplot() + 
  xlab('number of genomes') + 
  background_grid()
  # theme_bw()




sero %>% 
  select(asm_acc, Serotype) %>% 
  left_join(meta) %>%
  mutate(serotype=fct_lump((Serotype), n = 10)) %>%
  group_by(serotype, STX) %>% 
  tally() %>% 
  arrange(desc(n)) %>%
  # mutate(serotype=fct_inorder(serotype)) %>%
  ggplot(aes(x=STX, y=n)) + 
  geom_col(aes(fill=serotype)) + 
  ggtitle('serotypes of stx containing E.coli') + 
  # theme_cowplot() + 
  theme_bw()


# meta$asm_acc
# contains all SNP clusters that contain some O157:H7
O157_H7 <- 
  meta %>%
  filter(asm_acc != 'NULL') %>% 
  filter(PDS_acc %in% O157_PDS) %>%
  write_tsv('./output/O157:H7_meta.tsv')

# O157_H7 %>% 
#   pull(asm_acc) %>%
#   paste0(.,'.fna') %>%
#   write_lines('./data/genomes/O157_genomes.txt')


move_tib <- 
  tibble(files=list.files('./data/genomes', pattern = 'fna', full.names = T), 
         asm_acc=sub('./data/genomes/(.*).fna','\\1',files), 
         O157=ifelse(asm_acc %in% O157_H7$asm_acc, T, F)) %>% 
  mutate(dest_path=
           case_when(
           O157  ~ paste0('./data/O157/', asm_acc, '.fna'), 
           !O157 ~ paste0('./data/non_O157/', asm_acc, '.fna')
         ))

move_res <- file.rename(from = move_tib$files, to = move_tib$dest_path)

all(move_res)

# build ppanggolin list file


# want to indicate all complete genome contigs are circular
# these are complete genomes...
complete_genomes <- O157_H7[O157_H7$asm_level == 'Complete Genome',]


grep_vector <- function(pattern_vector, search_vector, INVERT=F){
  
  matches <- unique(grep(paste(pattern_vector,collapse="|"), 
                          search_vector,invert = INVERT, value=TRUE))
  return(matches)
  
}

complete_genome_paths <- 
  grep_vector(list.files('./data/O157', full.names = T),
              pattern_vector = complete_genomes$asm_acc)

infection_genomes_paths <- grep(list.files(path="./data/O157", full.names = T), pattern='GCA_', invert=TRUE, value=TRUE)


complete_genome_paths <- c(complete_genome_paths, infection_genomes_paths)



incomplete_genome_paths <-  
  grep_vector(list.files('./data/O157', full.names = T),
              pattern_vector = complete_genome_paths, 
              INVERT = T)

# RAN RENAME CONTIGS HERE ###


build_ppanggolin_file_fastas(complete_genome_paths = complete_genome_paths,
                             incomplete_genome_paths = incomplete_genome_paths) %>% 
  write_tsv(col_names = FALSE, file = 'ppangg_file.tsv')








