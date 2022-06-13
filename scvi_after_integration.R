library(dplyr)
library(Seurat)
library(ggplot2)
#library(colorspace)
#library(STACAS)
library(clustree)
library(gridExtra)
library(ggpubr)
#library(egg)
#library(scales)
library(scran)
#library(effectsize)
library(reticulate)
library(sceasy)
library(SingleCellExperiment)
library(scater)

#scanorama <- import('scanorama') 

print("importing...")
sc <- import("scanpy")
print("imported scanpy")
scvi <- import("scvi")
print("imported scvi")


filepath <- '~/projects/def-gbader/cvolk/'

samples <- readRDS(file=paste(filepath,'outputs/scvi_samples.RDS',sep=""))

pdf(file=paste(filepath,'outputs/immune_plots_scvi.pdf',sep=""), width=15, height=9)


samples <- FindNeighbors(samples, reduction="scvi", dims=1:10)
samples <- FindClusters(samples, resolution=0.2)
samples <- RunUMAP(samples, dims=1:10, reduction="scvi", n.components=2)

merge.immune <- samples

print("before dimplot")

DimPlot(merge.immune,
        group.by='SampleID')
FeaturePlot(merge.immune, features=c('CD3D', 'CD4', 'CD8A', 'CD44', 'CD68', 'ITGAM', 'CD274', 'CD28', 'PDCD1'), order=T)
DimPlot(merge.immune,
        group.by='seurat_clusters', label=T)
DimPlot(merge.immune, group.by='orig.ident')

dev.off()

pdf(file=paste(filepath,'outputs/immune_plots_scvi_split.pdf',sep=""),width=21, height=2)

DimPlot(merge.immune, split.by='orig.ident')

dev.off()

saveRDS(merge.immune, file=paste(filepath,'outputs/scvi_immune.RDS',sep=""))
