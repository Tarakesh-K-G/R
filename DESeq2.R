setwd("C:\\Users\\user\\Desktop\\RNA sequencing\\DESeq2")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("airway")
library(airway)
library(tidyverse)
data(airway)
airway
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)



#Step-1: Preparing count data......
#Read in counts data
counts_data<-read.csv("counts_data.csv")
head(counts_data)

#Read in sample info
colData<-read.csv("sample_info.csv")

#Making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data)%in%rownames(colData))

#are they in same order
all(colnames(counts_data)==rownames(colData))

#Step-2: Construct DESeq2 dataset object.......

dds<-DESeqDataSetFromMatrix(countData=counts_data,
                       colData=colData,
                       design = ~dexamethasone)
dds


#Pre-filtering: removing rows with low gene counts
#Keeping rows that have at least 10 reads total

keep<-rowSums(counts(dds))>=10
dds<- dds[keep,]
dds


#Set the factor level

dds$dexamethasone<-relevel(dds$dexamethasone,ref="untreated")

#NOte: Collapse technical replicates

#Run DESeq2.......
dds<-DESeq(dds)
res<-results(dds)
res

#Explore Result.......
summary(res)

#Adjust the p-value....
res0.01<-results(dds,alpha=0.01)
summary(res0.01)

#Contrast
resultsNames(dds)

#e.g.: treated_4hrs, treated_8hrs, untreated
results(dds,contrast=c("dexamethasone","treated_4hrs","untreated"))

#MA plot
plotMA(res)
