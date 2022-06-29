library(readxl)
library(pdtools)
library(tidyverse)
library(lubridate)

grep_vector <- function(pattern_vector, search_vector, INVERT=F){
  matches <- unique(grep(paste(pattern_vector,collapse="|"), 
                         search_vector,invert = INVERT, value=TRUE))
  return(matches)
  
}


usethis::use_directory('Nathan_requests')
# requesting isolates to phenotype and evaluate #

meta <- 
  read_tsv('output/O157:H7_meta.tsv') %>%
  filter(continent == 'Americas')

meta$collection_agency %>% unique()


NARMS_ISOLATES <- 
  read_xlsx('O157_NARMS.xlsx') %>%
  pull(`NCBI Accession Number`) %>%
  unique()

maybe_requestable <- 
  meta %>%
  filter(biosample_acc %in% NARMS_ISOLATES |   # NARMS freezer isolates 
           grepl('FSIS', collection_agency) |  # FSIS isolates 
           grepl('USDA', collection_agency) |  # other USDA isolates 
           grepl('FDA', collection_agency) |   # FDA isolates 
           State != '') %>%                    # Isolates with state level data 
  filter(Year >= 2018) %>%                     # from the last 5 years
  unique() %>%
  mutate(NARMS=ifelse(biosample_acc %in% NARMS_ISOLATES, TRUE, FALSE)) %>% 
  write_tsv('Nathan_requests/maybe_requestable_metadata.tsv')

maybe_requestable %>%
  pull(collected_by) %>%
  table()

set.seed(7)
maybe_requestable_summary_tree <- 
  maybe_requestable %>% 
  filter(!is.na(PDS_acc)) %>%
  group_by(PDS_acc) %>% 
  summarise(num_genomes=n(), 
            representative=asm_acc[sample(1:n(), 1)]) %>%
  arrange(desc(num_genomes)) %>% 
  write_tsv(paste0('Nathan_requests/maybe_requestable_tree_reps_',today(),'.tsv'))


# Construct ppanggolin file for requestable isolate representatives #

incomplete_genomes <- grep_vector(pattern_vector = maybe_requestable_summary_tree$representative, list.files('./data/O157', full.names = T))
complete_genomes <- grep('GCA_', list.files('./data/O157', full.names = T), invert = T, value = T)


build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes, 
                             incomplete_genome_paths = incomplete_genomes) %>% 
  write_tsv(paste0('maybe_requestable_reps',today(),'_ppanggolin.tsv'), col_names = FALSE)



######### copy all rep genomes to folder for ectyping ######

usethis::use_directory('Nathan_requests/dan_lin/')

dan_lin_copy <- 
  tibble(PATH=c(complete_genomes, incomplete_genomes)) %>% 
  mutate(NEW_PATH=sub('./data/O157/(.*)','./Nathan_requests/dan_lin/\\1',PATH))


if(!all(file.exists(dan_lin_copy$NEW_PATH))){
  
  copied <- file.copy(from = dan_lin_copy$PATH, to = dan_lin_copy$NEW_PATH, overwrite = TRUE)
  all(copied)
  
}



##### now one big pangenome for gene presence absence matrix ###
# LOOK <- maybe_requestable %>% filter(asm_acc == 'GCA_020997995.1')

# weird duplication going on...

maybe_requestable %>% filter(asm_acc == 'GCA_020757275.1') 


incomplete_genomes <- paste0('./data/O157/', maybe_requestable$asm_acc, '.fna')
complete_genomes <- grep('GCA_', list.files('./data/O157', full.names = T), invert = T, value = T)


build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes, 
                             incomplete_genome_paths = incomplete_genomes) %>% 
  write_tsv(paste0('maybe_requestable_all',today(),'_ppanggolin.tsv'), col_names = FALSE)


### now all the genomes....
usethis::use_directory('Nathan_requests_all/dan_lin/')

dan_lin_copy <- 
  tibble(PATH=c(complete_genomes, incomplete_genomes)) %>% 
  mutate(NEW_PATH=sub('./data/O157/(.*)','./Nathan_requests_all/dan_lin/\\1',PATH))


if(!all(file.exists(dan_lin_copy$NEW_PATH))){
  
  copied <- file.copy(from = dan_lin_copy$PATH, to = dan_lin_copy$NEW_PATH, overwrite = TRUE)
  all(copied)
  
}

# now run bash script


