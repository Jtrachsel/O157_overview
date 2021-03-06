

library(tidyverse)
library(cowplot)
library(Biostrings)

# ran ectyper on all downloaded genomes

# output of ectyper
sero <- read_tsv('./data/genomes/serotypes/output.tsv')
sero2 <- read_tsv('./data/genomes/2nd_serotypes/output.tsv')
sero <- bind_rows(sero, sero2)

# output of download script
meta <- read_tsv('./data/North_America_stx_metadata.tsv')


sero %>%
  mutate(serotype=fct_lump((Serotype), n = 10)) %>%
  group_by(serotype) %>% 
  tally() %>% arrange(desc(n)) %>%
  mutate(serotype=fct_inorder(serotype)) %>%
  ggplot(aes(x=serotype, y=n)) + geom_col() + 
  ggtitle('serotypes of stx containing North American E.coli from 2010 onward') + 
  # theme_cowplot() + 
  theme_bw()

# meta$asm_acc

O157_H7 <- 
  sero %>%
  mutate(asm_acc = sub('(GCA_[0-9]+.[0-9]+)_.*','\\1',Name)) %>%
  select(asm_acc, Serotype) %>%
  left_join(meta) %>%
  filter(Serotype == 'O157:H7') %>%
  write_tsv('./data/O157:H7_meta.tsv')

O157_H7 %>% 
  pull(asm_acc) %>%
  paste0(.,'.fna') %>%
  write_lines('./data/genomes/O157_genomes.txt')
# now I need to keep only O157 genomes?
# move O157 genomes to their own folder
## MADE DIRECTORIES AND MOVED TO NEW FOLDER ##




# build ppanggolin list file


# want to indicate all complete genome contigs are circular
# these are complete genomes...
complete_genomes <- O157_H7[O157_H7$asm_level == 'Complete Genome',]

# complete_genomes$asm_acc











grep_vector <- function(pattern_vector, search_vector, INVERT=F){
  
  matches <- unique(grep(paste(pattern_vector,collapse="|"), 
                          search_vector,invert = INVERT, value=TRUE))
  return(matches)
  
}

build_ppanggolin_file_fastas <-
  function(complete_genomes_paths=NULL,
           incomplete_genome_paths=NULL){
    if (is.null(complete_genomes_paths) & is.null(incomplete_genome_paths)){

      errorCondition('you must specify either complete_genome_paths or incomplete_genome_paths')
    } 
    complete_genomes_file <- 
      tibble(paths=complete_genome_paths, 
             ID=sub('.?/?.?/(.*)\\.f.*a$','\\1',paths), 
             fasta=map(.x = paths, .f=readDNAStringSet), 
             c_names=map(.x=fasta, .f=names), 
             contig_names=map_chr(.x=c_names, .f=~paste(.x, collapse='\t'))) |> 
      select(ID, paths, contig_names)
    
  }

complete_genome_paths <- 
  grep_vector(list.files('./data/genomes/O157', full.names = T),
              pattern_vector = complete_genomes$asm_acc)

infection_genomes_paths <- grep(list.files(path="./data/genomes/O157", full.names = T), pattern='GCA_', invert=TRUE, value=TRUE)


complete_genome_paths <- c(complete_genome_paths, infection_genomes_paths)

# 
# complete_genomes_file <- 
#   tibble(paths=complete_genome_paths, 
#        ID=sub('./data/genomes/O157/(.*)\\.f.*a','\\1',paths), 
#        fasta=map(.x = paths, .f=readDNAStringSet), 
#        c_names=map(.x=fasta, .f=names), 
#        contig_names=map_chr(.x=c_names, .f=~paste(.x, collapse='\t'))) |> 
#   select(ID, paths, contig_names)
# 
# complete_genomes_file


incomplete_genome_paths <-  
  grep_vector(list.files('./data/genomes/O157', full.names = T),
              pattern_vector = complete_genome_paths, 
              INVERT = T)

# incomplete_genomes_file <- 
#   tibble(paths=incomplete_genome_paths, 
#          ID=sub('./data/genomes/O157/(.*).fna','\\1',paths)) |> 
#   select(ID, paths)
# 


# ppanggolin_file <- 
#   bind_rows(complete_genomes_file, incomplete_genomes_file) |> 
  


# RAN RENAME CONTIGS HERE ###


build_ppanggolin_file_fastas(complete_genome_paths = complete_genome_paths,
                             incomplete_genome_paths = incomplete_genome_paths) %>% 
  write_tsv(col_names = FALSE, file = 'ppangg_file.tsv')








