setwd("C:\\Users\\user\\Desktop\\Single cell")
library(Seurat)
library(tidyverse)
A<-Read10X_h5(filename = 'kidney.h5')
str(A)
A[1:10,1:10]
Seurat.obj<-CreateSeuratObject(counts = A,project = "Kidney gene expression", 
                               min.cells = 5,min.features = 200)
Seurat.obj
str(Seurat.obj)

#Step 1.QC.....
View(Seurat.obj@meta.data)
#%MT reads[mitochondrial gene]
Seurat.obj[["percent.mt"]]<-PercentageFeatureSet(Seurat.obj,pattern = "^MT_")
View(Seurat.obj@meta.data)
VlnPlot(Seurat.obj,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
FeatureScatter(Seurat.obj,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+geom_smooth(method = 'lm')

#Step 2.Filtering.........
Seurat.obj<-subset(Seurat.obj,subset = nFeature_RNA>200&nFeature_RNA<2500)

#Step 3.Normalize data.....
Seurat.obj<-NormalizeData(Seurat.obj)
str(Seurat.obj)

#Step 4.Identify highly variable features......
Seurat.obj<-FindVariableFeatures(Seurat.obj,selection.method = "vst",nfeatures = 2000)

#Identify the most 10 highly variable genes..
top10<-head(VariableFeatures(Seurat.obj),10)

#Plot variable feature with & without labels...
plot1<-VariableFeaturePlot(Seurat.obj)
LabelPoints(plot = plot1, points = top10,repel = TRUE)

#Step 5.Scaling.....
all.genes<-rownames(Seurat.obj)
Seurat.obj<-ScaleData(Seurat.obj,features = all.genes)
str(Seurat.obj)

#Step 6.Perform Linear dimensionality reduction.....
Seurat.obj<-RunPCA(Seurat.obj,features = VariableFeatures(object = Seurat.obj))

#Visualize the PCA Result
print(Seurat.obj[["pca"]],dims = 1:5,nfeatures = 5)
DimHeatmap(Seurat.obj,dims = 1,cells = 500,balanced = TRUE)
# Determine dimensionality of the data....
ElbowPlot(Seurat.obj)

#Step 7.Clustering...
Seurat.obj<-FindNeighbors(Seurat.obj,dims = 1:20)

#Understanding resolution...
Seurat.obj<-FindClusters(Seurat.obj,resolution = c(0.3,0.5,0.7,1))
View(Seurat.obj@meta.data)

DimPlot(Seurat.obj,group.by = "RNA_snn_res.0.3",label = TRUE)

#Setting identity of clusters..
Idents(Seurat.obj)<-"RNA_snn_res.0.3"
Idents(Seurat.obj)
#Non-linear dimensionality reduction...
# If you haven't installed UMAP, you can do so via
reticulate::py_install(packages = 'umap-learn')

Seurat.obj<-RunUMAP(Seurat.obj,dims = 1:15)

#Note that you can set 'label=TRUE' or use the Label Clusters function to help label.
#Individual Clusters
DimPlot(Seurat.obj,reduction = "umap")
