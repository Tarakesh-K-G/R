setwd("C:\\Users\\user\\Desktop\\Doublets Finder")
library(Seurat)
library(ggplot2)
library(tidyverse)
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)

#Creat count matrix
cts<-ReadMtx(mtx = 'C:\\Users\\user\\Desktop\\Doublets Finder\\raw_feature_bc_matrix\\matrix.mtx.gz',
             features = 'C:\\Users\\user\\Desktop\\Doublets Finder\\raw_feature_bc_matrix\\features.tsv.gz',
             cells = 'C:\\Users\\user\\Desktop\\Doublets Finder\\raw_feature_bc_matrix\\barcodes.tsv.gz')
cts[1:10,1:10]

#Creat Seurat Object
pbmc.seurat<-CreateSeuratObject(counts=cts)
str(pbmc.seurat)

#Qc and Filtering
#explore Qc
pbmc.seurat$mtpercent<-PercentageFeatureSet(pbmc.seurat,pattern = '^MT-')
#Filtering
pbmc.seurat.filterd<-subset(pbmc.seurat,subset = nCount_RNA>800 & nFeature_RNA>500 & mtpercent<10)
pbmc.seurat
pbmc.seurat.filterd

#pre-process workflow
pbmc.seurat.filterd<-NormalizeData(object = pbmc.seurat.filterd)
pbmc.seurat.filterd<-FindVariableFeatures(object = pbmc.seurat.filterd)
pbmc.seurat.filterd<-ScaleData(object = pbmc.seurat.filterd)
pbmc.seurat.filterd<-RunPCA(object = pbmc.seurat.filterd)
ElbowPlot(pbmc.seurat.filterd)
pbmc.seurat.filterd<-FindNeighbors(object = pbmc.seurat.filterd, dims = 1:20)
pbmc.seurat.filterd<-FindClusters(object = pbmc.seurat.filterd)
pbmc.seurat.filterd<-RunUMAP(object = pbmc.seurat.filterd,dims = 1:20)

#pK Identification (no ground-truth)----
sweep.res.list_pbmc <- paramSweep(pbmc.seurat.filterd, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc.seurat.filterd@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(pbmc.seurat.filterd@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
pbmc.seurat.filterd<-doubletFinder(pbmc.seurat.filterd, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")


# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)