library(pdtools)
library(tidyverse)
library(rvest)
library(lubridate)
library(curl)
library(data.table)
library(cowplot)

system('mkdir ./data/')

pdtools::list_PDGs('Escherichia_coli_Shigella')

# downloads the most recent complete NCBI path detection metadata
# download_most_recent_complete('Escherichia_coli_Shigella', folder_prefix = './data/')



# download_PDD_metadata('Escherichia_coli_Shigella', PDG = 'PDG000000004.2937', folder_prefix = './data/' )


# curl_download('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt', 
              # destfile = './data/assembly_summary.txt')

ass_sum <- fread('./data/assembly_summary.txt', quote = '') %>%
  mutate(asm_acc=`# assembly_accession`)

AMR <-  fread('data/PDG000000004.2937.amr.metadata.tsv', quote = '')
clusts <- fread('data/PDG000000004.2937.cluster_list.tsv', quote = '')

master_meta <- AMR %>% left_join(clusts) 

# detect broad STX types
master_meta <- 
  master_meta %>%
  mutate(STX=case_when(
  grepl('stxA', virulence_genotypes) & grepl('stxB', virulence_genotypes)  ~ 'both', 
  grepl('stxA', virulence_genotypes) ~ 'stxA', 
  grepl('stxB', virulence_genotypes) ~ 'stxB', 
  TRUE ~ 'none'
))


# within PDS accessions there is some diversity in the presence of STX genes

master_meta %>% 
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
master_meta %>% 
  select(target_acc, virulence_genotypes) %>%
  filter(grepl('stx', virulence_genotypes)) %>%
  mutate(vir = gsub('\"','',virulence_genotypes)) %>%
  select(-virulence_genotypes) %>%
  separate_rows(vir, sep = ',') %>%
  filter(grepl('stx', vir)) %>%
  group_by(vir) %>% tally() %>%
  filter(!grepl('=', vir)) %>%
  mutate(type=sub('(stx[A-B][1-2])[a-z]','\\1',vir)) %>%
  ggplot(aes(x=vir, y=n)) +geom_col()+ facet_wrap(~type, scales = 'free')

         
##
PDD_meta <- 
  master_meta %>%
  get_earliest_year() %>%
  mutate(country=sub('([A-Za-z ]+):([A-Za-z ]+)','\\1',geo_loc_name)) 

  
PDD_meta <- 
  PDD_meta %>% 
  filter(grepl('USA|Canada|Mexico', country, ignore.case = T)) %>%
  filter(epi_type != 'NULL') %>%
  filter( Year > 2010)
#####



genomes_per_year_North_America <- 
  PDD_meta %>%  
  count(Year, epi_type, name = 'tot_epi_yr')
  

  


p1 <- 
  PDD_meta %>%
  group_by(STX, epi_type) %>%
  tally()  %>%
  ggplot(aes(x=STX, y=n, fill=STX)) + 
  geom_col(position = position_dodge(), color='black', show.legend = F) +
  geom_text(aes(label=scales::comma(n, accuracy = 1), group=STX),color='black', nudge_y = 1000)+ 
  facet_wrap(~epi_type, ncol = 1, scales = 'free_x') +
  theme_cowplot()+
  theme(panel.grid.major = element_line(color='grey'), 
        panel.border = element_rect(color='black'), 
        axis.title = element_text(size=16)) +
  ylab('number of genomes') + 
  xlab('stx gene presence')



p2 <- 
  PDD_meta %>%
  group_by(STX, Year, epi_type) %>%
  tally() %>% 
  left_join(genomes_per_year_North_America) %>%
  mutate(perc_epi_yr = (n / tot_epi_yr)*100 ) %>% 
  filter(epi_type != 'NULL') %>% 
  ggplot(aes(x=Year, y=perc_epi_yr, color=STX)) + 
  geom_line(size=1)+
  geom_point(aes(fill=STX), color='white',size=3, shape=21, stroke=1.25) +
  facet_wrap(~ epi_type, ncol = 1, scales = 'free_x') + 
  theme_cowplot() + 
  theme(panel.grid.major = element_line(color='grey'), 
        panel.border = element_rect(color='black'), 
        axis.title = element_text(size=16)) + 
  ylab('percent isolates per year') #+ 


# p1
# p2

plot_grid(p1, p2, 
          labels = c('A', 'B'),
          rel_widths = c(1,2) ,
          label_size = 16,
          nrow = 1)

stx_meta <- PDD_meta %>% filter(STX == 'both')


stx_meta %>%
  group_by(Year, epi_type) %>%
  tally() %>%
  ggplot(aes(x=Year, y=n, color=epi_type)) +
  geom_point() +
  geom_line() + 
  theme_cowplot()+
  theme(panel.grid= element_line(color='grey'), 
        panel.border = element_rect(color='black'), 
        axis.title = element_text(size=16)) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  ggtitle('genomes containing stxA and stxB')


stx_meta %>% write_tsv('./data/North_America_stx_metadata.tsv')


PDD_meta %>%
  filter(STX == 'both') %>%
  filter(grepl('USA|Canada|Mexico', country, ignore.case = T)) %>%
  count(serovar) %>% 
  arrange(desc(n))


O157_PDSs <- 
  PDD_meta %>%
  filter(grepl('O157', serovar)) %>%
  pull(PDS_acc) %>%
  unique()

stx_meta %>% 
  filter(PDS_acc %in% O157_PDSs) %>%
  filter(STX == 'both')



ass_sum <- ass_sum %>% filter(asm_acc %in% stx_meta$asm_acc)

# ass_sum$ftp_path

system('mkdir ./data/genomes/')



make_fna_urls(stx_meta$asm_acc, assembly_summary = ass_sum) %>% 
  write_lines('./data/genomes/North_America_stx.ftp_paths')


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
