library(tidyverse)
library(dplyr)

setwd("C:\\Users\\user\\Desktop\\RNA sequencing")
Data<-read.csv('GSE183947_fpkm.csv')
dim(Data)

#Get metadata
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GEOquery")
library(GEOquery)
gse<-getGEO(GEO = 'GSE183947',GSEMatrix=TRUE)
gse
metadata<-pData(phenoData(gse[[1]]))
metadata
head(metadata)
View(metadata)

metadata.modified<-metadata %>% select(c(1,10,11,17)) %>% 
  rename(tissue=characteristics_ch1,metastasis=characteristics_ch1.1) %>% 
mutate(tissue=gsub("tissue:","",tissue),metastasis=gsub("metastasis:","",metastasis))
metadata.modified
#gsub:global substitution
head(Data)
#reshaping Data
Data.long<-Data %>% rename(gene=X) %>% gather(samples,FPKM,-gene)

#Join dataframes = Data.long + metadata.modified

ABC.long<-Data.long %>% left_join(.,metadata.modified,by=c('samples'='description'))
ABC.long

#explore data

ABC.long %>% filter(gene=='BRCA1'|gene=='BRCA2') %>% 
  group_by(gene,tissue) %>% summarize(mean_FPKM=mean(FPKM)) %>% arrange(mean_FPKM)
