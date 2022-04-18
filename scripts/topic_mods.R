#Think about novel gene pairs
# how often do pairs of gene co-occur?


library(pdtools)
library(tidyverse)

# pan_PA <- generate_pangenome(num_genomes = 100, num_genes = 7000, core_genome_fraction = .7)


# these co-occurence matricies get really large really fast...
# need some heuristics to reduce the problem?

# remove persistent partition?
library(data.table)

pan_PA <- fread('../O157_overview/pan/gene_presence_absence.Rtab', quote='')
shell <- read_lines('../O157_overview/pan/partitions/shell.txt')
cloud <- read_lines('../O157_overview/pan/partitions/cloud.txt')

pan_PA <-
  pan_PA %>%
  filter(Gene %in% cloud) %>%
  column_to_rownames(var = 'Gene') %>%
  as.matrix() %>% 
  as('sparseMatrix')


summary(pan_PA)

pan_PA <- pan_PA[rowSums(pan_PA)/nrow(pan_PA) > .005,]

# I think this removes columns with all zeros?
pan_PA <- pan_PA[,unique(summary(pan_PA)$j)] # returns non_zero column indexes, then removes



library(Matrix)
M <- Matrix(c(0,0,0,1,0,0,0,1,1,1,0,0), nc=4)
M[,unique(summary(M)$j)]





# pan_PA[,sample(1:ncol(pan_PA))]

library(Matrix)

CP <- Matrix::tcrossprod(pan_PA)
diag(CP) <- 0
CP <- drop0(CP)

# dim(CP)
# dim(CP_sparse_orig)

# sum(CP_sparse_orig == 0)/length(CP_sparse_orig)
sum(CP == 0)/length(CP)

one_genome_gene_adj <- tcrossprod(pan_PA[,1, drop = FALSE])

test <- apply(CP[,1:3], 2, Matrix::tcrossprod, simplify = FALSE)
whatsthis <- test[[1]]

m <- pan_PA

listCols<-function(m){
  # browser()
  #converts a sparse Matrix into a list of its columns
  res<-split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
  names(res)<-colnames(m)
  res
}


listCols2<-function(m){
  #using slam simple sparse array
  m2<-as.simple_sparse_array(m)
  res<-split(m2$v,m2$i[,2])
  names(res)<-colnames(m)
  res
}

listCols3<-function(m){
  #using slam simple sparse triplet
  m2<-as.simple_triplet_matrix(m)
  res<-split(m2$v,m2$j)
  names(res)<-colnames(m)
  res
}

listCols4<-function(m){
  #using data.frame
  m2<-summary(m)
  res<-split(m2$x,m2$j)
  names(res)<-colnames(m)
  res
}

listCols5<-function(m){
  #using data.table
  m2<-as.data.table(summary(m))
  res<-split(m2$x,m2$j)
  names(res)<-colnames(m)
  res
}

res<-microbenchmark(listCols(m),
                    listCols2(m),
                    listCols3(m),
                    listCols4(m),
                    listCols5(m),
                    unit="s",times=10)
boxplot(res,unit="s")

tst <- listCols(pan_PA)
names(tst[[1]])

baseline<-function(m,f,nproc=2){
  unlist(mclapply(listCols(m),f,mc.cores=nproc))
}


map


CP[1:5, 1:5]
min(CP)
MAX_COC <- max(CP)
hist(log(CP))

# looking to select genomes with many novel gene combinations
# so low co-occurence gene_pairs are of interest
# give each genome a score based on novel gene combos?
  # num gene combos below a certain threshold?
library(tidyverse)

CP_long <-
  CP %>%
  as.data.frame() %>%
  rownames_to_column(var = 'GENE1') %>%
  pivot_longer(cols = -GENE1, names_to = 'GENE2', values_to = 'num_co_occur') %>%
  filter(num_co_occur > 0) %>%
  mutate(prop_something=num_co_occur/MAX_COC)

CP[1:2,1:2]




hist(CP_long$prop_something)
CP_long %>% arrange(desc(num_co_occur))
library(fastTopics)

library(Matrix)
meta <- read_tsv('../O157_overview/data/O157:H7_meta.tsv')
pan_PA <- as(pan_PA, "sparseMatrix")


# each genome is a mixture of topics
# each gene is a mixture of topics
topmod <- fit_topic_model(pan_PA, k=4, control.main = list(numiter = 4,extrapolate = TRUE,nc = 40))

# topics per genome
meta <- topmod$F %>% as.data.frame %>% rownames_to_column(var = 'asm_acc') %>% right_join(meta)

meta %>% select(ag_match)

meta %>% select(Year, starts_with('k')) %>%
  pivot_longer(-Year, names_to = 'topic', values_to = 'prop') %>%
  group_by(Year, topic) %>%
  summarise(topic_per_year=prop/sum(prop, na.rm = T)) %>%
  ggplot(aes(x=Year, y=topic_per_year, color=topic)) + geom_point()


meta$Year %>% max()
meta_topics <- meta %>% select(Year,asm_acc, starts_with('k')) %>%
  pivot_longer(-c(Year, asm_acc), names_to = 'topic', values_to = 'prop')


meta_topics[which(is.na(meta_topics$prop)),]

topmod$F[rownames(topmod$F) == 'GCA_011818065.1',]
# topics per gene
topmod$L


rowSums(topmod$L)
colSums(topmod$L)

topmod$Fn
topmod$Fy

