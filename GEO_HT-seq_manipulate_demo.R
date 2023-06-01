# this script is to manipulate high throughout RNA-seq data downloaded from GEO
# GSE183947

# load packages------
library(dplyr)
library(tidyverse)
library(GEOquery)

# read in expression data------
# set direction
dir <- "/Users/wanzhihan/Desktop/R语言/05_GEO转录组视频课程代码/demo/"
dat <- read.csv(file = paste(dir,"02_data input/GSE183947_fpkm.csv", sep=""))

# check the dimension of your dataframe
dim(dat) # 20246    61

# get metadata using geoquery------
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)

# if the connection buffer was not large enough
# run: Sys.setenv("VROOM_CONNECTION_SIZE" = ......*1000)

# quickly check the content
gse

# get the phynoData from gse list 1(metadata)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get interested information
metadata.subset <- select(metadata, c(1,10,11,17))
head(metadata.subset)
dim(metadata.subset)
# lets do it using dplyr
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue:","", tissue)) %>%
  mutate(metastasis = gsub("metastasis:", "", metastasis))

#add those information into gene expression data------
head(dat)

# reshaping data
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene) 
  
# join dataframes = dat.long + metadata.modified
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

#examples of exploring data------
dat.long %>%
  filter(gene == 'BRCA1'|gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM), median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)
