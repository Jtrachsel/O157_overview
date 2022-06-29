library(pdtools)
library(tidyverse)
library(data.table)
library(cowplot)
library(usethis)
library(scales)
library(countrycode)


# system('mkdir ./data/')
use_directory('./data/')
use_directory('./data/genomes')
use_directory('./data/O157')
use_directory('./data/non_O157')
# use_directory('./data/genomes/serotypes/')

pdtools::list_PDGs('Escherichia_coli_Shigella')

# downloads the most recent complete NCBI path detection metadata

#TODO add check if most recent complete exists already
dlmrc <- 
  download_most_recent_complete('Escherichia_coli_Shigella',
                                 folder_prefix = './data/')

download_gbk_assembly_summary('./data/assembly_summary.txt')


##### Read in most recent version of the metadata
# remove old metadata


# wrap this in function
# designed to clean up old versions of PDG metadata
# get_PDG_version('./data/')

most_recent_PDG <- 
  tibble(files=list.files('./data/', 'PDG.*.amr.metadata.tsv'), 
       PDG_num=as.numeric(sub('PDG([0-9]+.[0-9]+).amr.metadata.tsv','\\1',files)), 
       PDG=sub('(PDG[0-9]+.[0-9]+).amr.metadata.tsv','\\1',files)) %>% 
  arrange(desc(PDG_num)) %>% 
  slice_head(n=1) %>% 
  pull(PDG)

OLD_PDG_files <- grep(list.files(path="data", 'PDG', full.names = T), pattern=most_recent_PDG, invert=TRUE, value=TRUE)

most_recent_PDG_files <- list.files(path = './data/',pattern =  most_recent_PDG, full.names = T)

removed <- file.remove(OLD_PDG_files) 

if (!all(removed)){
  print(paste('These files could be removed:'))
  print(OLD_PDG_files[!removed])
}

###
meta <- 
  read_tsv(most_recent_PDG_files[1]) %>% 
  left_join(read_tsv(most_recent_PDG_files[2]))

# 

ass_sum <- fread('./data/assembly_summary.txt', quote = '') %>%
  mutate(asm_acc=`# assembly_accession`)

## check for genomes that have not had 'virulence' genes detected
# any(is.na(meta$virulence_genotypes))
# any(meta$virulence_genotypes == 'NULL')
# any(is.null(meta$virulence_genotypes))

# remove genomes that haven't had virulence genes detected
meta <- 
  meta %>% 
  filter(virulence_genotypes != 'NULL') %>% 
  filter(epi_type != 'NULL')


### classify each genome in terms of a broad STX type
meta <- 
  meta %>%
  mutate(
    STX=case_when(
    grepl('stxA1', virulence_genotypes) &
      grepl('stxB1', virulence_genotypes) &
      grepl('stxA2', virulence_genotypes) &
      grepl('stxB2', virulence_genotypes)  ~ 'Both',
    grepl('stxA1', virulence_genotypes) &
      grepl('stxB1', virulence_genotypes)  ~ 'STX1', 
    grepl('stxA2', virulence_genotypes) &
      grepl('stxB2', virulence_genotypes)  ~ 'STX2', 
    grepl('stxA1', virulence_genotypes) ~ 'one subunit',
    grepl('stxB1', virulence_genotypes) ~ 'one subunit',
    grepl('stxA2', virulence_genotypes) ~ 'one subunit',
    grepl('stxB2', virulence_genotypes) ~ 'one subunit',
    TRUE ~ 'None'
  ),
    STX=factor(STX, levels = c('None', 'Both', 'STX1', 'STX2')))


# within PDS accessions there is some diversity in the presence of STX genes
# This only looks at the presence of both subunits not stx types
meta %>% 
  group_by(PDS_acc, STX) %>%
  tally() %>% 
  tally() %>%
  arrange(desc(n)) %>% 
  count(n) %>%
  ggplot(aes(x=n, y=nn)) + geom_col() + 
  ggtitle('number of different STX types in PDS groups', 
          'most PDS groups only have one STX type')



### what are all the stx genes?

# stx subtypes
meta %>% 
  select(target_acc, virulence_genotypes) %>%
  filter(grepl('stx', virulence_genotypes)) %>%
  mutate(vir = gsub('\"','',virulence_genotypes)) %>%
  select(-virulence_genotypes) %>%
  separate_rows(vir, sep = ',') %>%
  filter(grepl('stx', vir)) %>%
  group_by(vir) %>% tally() %>%
  filter(!grepl('=', vir)) %>%
  mutate(type=sub('(stx[A-B][1-2])[a-z]','\\1',vir)) %>%
  ggplot(aes(x=vir, y=n)) +
  geom_col() + 
  facet_wrap(~type, scales = 'free')+ 
  theme_half_open() + 
  background_grid() + 
  theme(axis.text.x=element_text(angle = -45, hjust = .1)) + 
  ylab('number of genomes')+
  xlab('gene')


ggsave('./output/stx_genes.jpeg', width = 8, height=5, bg='white')
         
###
# stx type for each genome

stx_types <- 
  meta %>% 
  select(target_acc, virulence_genotypes, STX) %>%
  filter(grepl('stx', virulence_genotypes)) %>%
  mutate(vir = gsub('\"','',virulence_genotypes)) %>%
  select(-virulence_genotypes) %>%
  separate_rows(vir, sep = ',') %>%
  filter(grepl('stx', vir)) %>% 
  filter(!grepl('=PARTIAL', vir)) %>% 
  filter(!grepl('=MISTRANSLATION', vir)) %>% 
  # mutate(vir=ifelse(grepl('=HMM', vir), 'stxHMM', vir)) %>% 
  filter(!grepl('=HMM', vir)) %>% 
  group_by(target_acc, STX) %>% 
  summarise(stx_type=paste(vir, collapse = '_'))


look <- 
  stx_types %>% 
  ungroup() %>% 
  mutate(type_lump=fct_lump_n(f = stx_type, n = 15))




cbp1 <- c(None="#999999", Both="#E69F00", STX1="#56B4E9", STX2="#009E73")

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# 
# 
# cbp1 <- c(None="#999999", Both="#E69F00", STX1="#56B4E9", STX2="#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


look %>%
  group_by(type_lump, STX) %>%
  tally() %>%
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(type_lump=fct_inorder(type_lump)) %>% 
  filter(!is.na(STX)) %>% 
  filter(type_lump != 'stxB1a') %>% 
  #filter(type_lump != 'Other') %>% 
  ggplot(aes(y=fct_rev(type_lump), x= n)) +
  geom_col(aes(fill=STX), color='black') + 
  theme_half_open() + 
  background_grid() + 
  scale_fill_manual(values = cbp1) + 
  ylab('STX type') + 
  xlab('number of genomes')
  

ggsave('./output/STX_subunit_types.jpeg', width = 8, height = 5, units = 'in', bg='white')

###

#### All worldwide isolates

year_dat <- meta %>% extract_earliest_year()
meta <- meta %>% left_join(year_dat)

genomes_per_year <-
  meta %>%
  count(Year, epi_type, name = 'tot_epi_yr')


# color blind palettes from:
# https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/



p1 <- 
  meta %>%
  filter(STX != 'one subunit') %>% 
  group_by(STX, epi_type) %>%
  tally()  %>%
  ggplot(aes(x=STX, y=n, fill=STX)) + 
  geom_col(position = position_dodge(), color='black', show.legend = F) +
  geom_text(aes(label=scales::comma(n, accuracy = 1), group=STX),color='black', nudge_y = 3500)+ 
  facet_wrap(~epi_type, ncol = 1, scales = 'free_x') +
  theme_cowplot()+
  scale_fill_manual(values = cbp1)+
  theme(panel.grid.major = element_line(color='grey'), 
        panel.border = element_rect(color='black'), 
        axis.title = element_text(size=16)) +
  ylab('number of genomes') + 
  xlab('stx gene presence')

p1


p2 <- 
  meta %>%
  filter(STX != 'one subunit') %>% 
  filter(Year > 2011) %>% 
  group_by(STX, Year, epi_type) %>%
  tally() %>% 
  left_join(genomes_per_year) %>%
  mutate(perc_epi_yr = (n / tot_epi_yr)*100 ) %>% 
  filter(epi_type != 'NULL') %>% 
  ggplot(aes(x=Year, y=perc_epi_yr, color=STX)) + 
  geom_line(size=1, )+
  geom_point(aes(fill=STX, shape=STX), color='white',size=3,stroke=1.25) +
  facet_wrap(~ epi_type, ncol = 1, scales = 'free_x') + 
  theme_cowplot() + 
  scale_shape_manual(values=c(None=21,Both=22, STX1=25, STX2=24))+
  scale_fill_manual(values = cbp1)+
  scale_color_manual(values = cbp1)+
  theme(panel.grid.major = element_line(color='grey'), 
        panel.border = element_rect(color='black'), 
        axis.title = element_text(size=16)) + 
  scale_x_continuous(breaks= pretty_breaks())+
  scale_y_continuous(breaks= pretty_breaks())+
  
  ylab('percent isolates per year') #+ 




p2

plot_grid(p1, p2, 
          labels = c('A', 'B'),
          rel_widths = c(1,2) ,
          label_size = 16,
          nrow = 1) + 
  ggtitle('stx toxin presence across all E. coli worldwide')

ggsave('./output/STX_over_time.jpeg', width = 9, height = 6, units = 'in', bg='white')

# we want to consider any genome with any combination of STX genes
stx_meta <- meta %>% filter(STX %in% c('Both', 'STX1', 'STX2'))

### a little metadata cleanup...

country_dat <- stx_meta %>% extract_country()
host_dat <- stx_meta %>% extract_consensus_ag_species()
agency_dat <- stx_meta %>% extract_collection_agency()


stx_meta <- 
  stx_meta %>% 
  # left_join(year_dat) %>% # commented out because I already joined in earlier
  left_join(country_dat) %>% 
  left_join(host_dat) %>% 
  left_join(agency_dat)



stx_meta <-
  stx_meta %>% 
  mutate(continent=countrycode(sourcevar = country, 
                               origin = 'country.name', 
                               destination = 'continent'))
## ADD extract continent

# stx_meta %>%
#   group_by(Year, epi_type, STX) %>%
#   tally() %>%
#   ggplot(aes(x=Year, y=n, color=epi_type)) +
#   geom_point() +
#   geom_line() + 
#   theme_cowplot()+
#   theme(panel.grid= element_line(color='grey'), 
#         panel.border = element_rect(color='black'), 
#         axis.title = element_text(size=16)) + 
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
#   ggtitle('genomes containing stx1, stx2, or both')


stx_meta %>% count(asm_acc == 'NULL')

stx_meta <- stx_meta %>% filter(asm_acc != 'NULL')



stx_meta %>% write_tsv('./output/all_stx_metadata.tsv')

# 
# meta %>%
#   # filter(STX == 'both') %>%
#   filter(grepl('USA|Canada|Mexico', country, ignore.case = T)) %>%
#   count(serovar) %>% 
#   arrange(desc(n))

# 
# O157_PDSs <- 
#   meta %>%
#   filter(grepl('O157', serovar)) %>%
#   pull(PDS_acc) %>%
#   unique()
# # 
# O157_meta <- stx_meta %>% 
#   filter(PDS_acc %in% O157_PDSs) 



# ass_sum <- ass_sum %>% filter(asm_acc %in% stx_meta$asm_acc)

# ass_sum$ftp_path

# system('mkdir ./data/genomes/')
# I downloaded all 43000+ genomes and used ectyper on them

# stx_meta <- stx_meta %>% filter(asm_acc != 'NULL')

stx_meta <- 
  stx_meta %>%
  make_ftp_paths('./data/assembly_summary.txt') %>% 
  make_download_urls(type = 'fna') %>% 
  make_dest_paths(type = 'fna', dest_dir = 'data/genomes/') #%>% 
  # download_genomes()

existing_assems <- 
  c(list.files('data/O157', 'GCA'), list.files('data/non_O157/', 'GCA'))
# list.files('data', 'GCA', recursive = T)
ex_asm_accs <- sub('.fna','',existing_assems)

download_these <-
  stx_meta %>%
  filter(!(asm_acc %in% ex_asm_accs)) %>% 
  download_genomes('fna')




# existing_asm_accs

# make_fna_urls(stx_meta$asm_acc, assembly_summary = ass_sum) %>% 
#   write_lines('./data/genomes/North_America_stx.ftp_paths')


# cat North_America_stx.ftp_paths | parallel -j 2 'wget {}'
###

###
# We compared the HUS risk between isolates with stx2a and those with stx2a and
# another gene and estimated additive interaction of the stx genes. Adjusted for
# age and symptoms, the HUS incidence of E. coli O157:H7 containing stx2a alone
# was 4.4% greater (95% confidence interval (CI) −0.3%, 9.1%) than when it
# occurred with stx1a. When stx1a and stx2a occur together, the risk of HUS was 
# 27.1% lower (95% CI −87.8%, −2.3%) than would be expected if interaction were
# not present
