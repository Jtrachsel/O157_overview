library(pdtools)
library(tidyverse)
library(data.table)
library(cowplot)
library(usethis)
library(scales)

# system('mkdir ./data/')
use_directory('./data/')


pdtools::list_PDGs('Escherichia_coli_Shigella')

# downloads the most recent complete NCBI path detection metadata
dlmrc <- download_most_recent_complete('Escherichia_coli_Shigella', folder_prefix = './data/')
download_gbk_assembly_summary('./data/assembly_summary.txt')


##### 

# wrap this in function
# designed to clean up old versions of PDG metadata


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


### Update collection function here

# generate_new_ftp_urls <- function(data_dir, organism, old_accessions){
#   dlmrc <- download_most_recent_complete(organism, folder_prefix = data_dir)
#   
#   
#   # 
#   tibble(meta_files=list.files(data_dir, pattern = 'PDG[0-9]+.[0-9]+.amr.metadata.tsv'), 
#          cluster_files=list.files(data_dir, pattern = 'PDG[0-9]+.[0-9]+.cluster_list.tsv')) %>% 
#     mutate(PDG=as.numeric(sub('PDG([0-9]+.[0-9]+).amr.metadata.tsv','\\1',meta_files))) %>% 
#     arrange(desc(PDG)) %>% 
#     select()
#   
#   
# }


# 
# 
# meta_files <- list.files('./data/', pattern = 'PDG[0-9]+.[0-9]+.amr.metadata.tsv')
# cluster_files <- list.files('./data/', pattern = 'PDG[0-9]+.[0-9]+.cluster_list.tsv')

#

# download_PDD_metadata('Escherichia_coli_Shigella', PDG = 'PDG000000004.2937', folder_prefix = './data/' )


# curl_download('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt', 
              # destfile = './data/assembly_summary.txt')

ass_sum <- fread('./data/assembly_summary.txt', quote = '') %>%
  mutate(asm_acc=`# assembly_accession`)

# AMR <-  fread('data/PDG000000004.2937.amr.metadata.tsv', quote = '')
# clusts <- fread('data/PDG000000004.2937.cluster_list.tsv', quote = '')

# meta <- AMR %>% left_join(clusts) 

# detect broad STX types
# meta <- 
#   meta %>%
#   mutate(STX=case_when(
#   grepl('stxA', virulence_genotypes) & grepl('stxB', virulence_genotypes)  ~ 'both', 
#   grepl('stxA', virulence_genotypes) ~ 'stxA', 
#   grepl('stxB', virulence_genotypes) ~ 'stxB', 
#   TRUE ~ 'none'
# ))

any(is.na(meta$virulence_genotypes))
any(meta$virulence_genotypes == 'NULL')
any(is.null(meta$virulence_genotypes))

# remove genomes that havent had virulence genes detected
meta <- meta %>% 
  filter(virulence_genotypes != 'NULL') %>% 
  filter(epi_type != 'NULL')



### The above makes sure the subunits are both tehre, but doesnt actually 
# differentiate the different types of stx, ie 1 and 2

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
  ggplot(aes(x=vir, y=n)) + geom_col() + facet_wrap(~type, scales = 'free')

         
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


look %>%
  group_by(type_lump, STX) %>% tally() %>%
  arrange((n)) %>% 
  ungroup() %>% 
  mutate(type_lump=fct_infreq(type_lump)) %>% 
  #filter(type_lump != 'Other') %>% 
  ggplot(aes(y=type_lump, x= n)) +
  geom_col(aes(fill=STX))



###

year_dat <- meta %>% extract_earliest_year()
country_dat <- meta %>% extract_country()
host_dat <- meta %>% extract_consensus_ag_species()
agency_dat <- meta %>% extract_collection_agency()
  

meta <- meta %>% 
  left_join(year_dat) %>% 
  left_join(country_dat) %>% 
  left_join(host_dat) %>% 
  left_join(agency_dat)

#### All worldwide isolates

# genomes_per_year_North_America <- 
#   meta %>%  
#   count(Year, epi_type, name = 'tot_epi_yr')
  

# color blind palettes from:
# https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
  left_join(genomes_per_year_North_America) %>%
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




stx_meta <- meta %>% filter(STX %in% c('Both', 'STX1', 'STX2'))


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


stx_meta %>% write_tsv('./data/North_America_stx_metadata.tsv')


meta %>%
  filter(STX == 'both') %>%
  filter(grepl('USA|Canada|Mexico', country, ignore.case = T)) %>%
  count(serovar) %>% 
  arrange(desc(n))


O157_PDSs <- 
  meta %>%
  filter(grepl('O157', serovar)) %>%
  pull(PDS_acc) %>%
  unique()

O157_metastx_meta %>% 
  filter(PDS_acc %in% O157_PDSs) %>%
  filter(STX == 'both')



# ass_sum <- ass_sum %>% filter(asm_acc %in% stx_meta$asm_acc)

# ass_sum$ftp_path

# system('mkdir ./data/genomes/')
# I downloaded all 43000+ genomes and used ectyper on them
stx_meta <- 
  stx_meta %>%
  make_ftp_paths('./data/assembly_summary.txt') %>% 
  make_download_urls(type = 'fna') %>% 
  make_dest_paths(type = 'fna', dest_dir = 'data/genomes/') #%>% 
  # download_genomes()

existing_assems <- 
  c(list.files('data/genomes/O157', 'GCA'), list.files('data/genomes/non_O157/', 'GCA'))
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
