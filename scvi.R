
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

m.immune <- readRDS(file=paste(filepath,'outputs/m_immune.RDS',sep=""))

### Initial cluster visualization ###

pdf(file=paste(filepath,'outputs/immune_plots_scvi.pdf',sep=""), width=15, height=9)

DimPlot(m.immune,
        group.by='SampleID')
FeaturePlot(m.immune, features=c('CD3D', 'CD4', 'CD8A', 'CD44', 'CD68', 'ITGAM', 'CD274', 'CD28', 'PDCD1'), order=T)
DimPlot(m.immune,
        group.by='seurat_clusters', label=T)

#DimPlot(m.tumor,
#        group.by='SampleID')
#FeaturePlot(m.tumor, features=c('PTPRC', 'CD24', 'ITGAM', 'HLA-A'), order=T)
#DimPlot(m.tumor,
#        group.by='seurat_clusters', label=T)


### Re-normalize data ###

#m.immune <- SCTransform(m.immune)

### Object -> list

sample.list <- SplitObject(m.immune, split.by = "orig.ident")
#sce.list <- lapply(sample.list, as.SingleCellExperiment)
#sce.list <- sample.list
### Re-normalization & Integration ###

#sce.list <- sample.list

for (i in 1:length(sample.list)) {
#  clusters <- quickCluster(sce.list[[i]], min.size=40)
#  sce.list[[i]] <- computeSumFactors(sce.list[[i]], clusters=clusters)
#  summary(sizeFactors(sce.list[[i]]))
#  sce.list[[i]] <- logNormCounts(sce.list[[i]])
#  print("norm counts")
 # DefaultAssay(sce.list[[i]]) <- "SCT"
 # print("default SCT")
  sample.list[[i]] <- SCTransform(sample.list[[i]])
#  print("SCTransform")
  sample.list[[i]] <- RunPCA(sample.list[[i]], npcs=30)
#  print("pca")
  sample.list[[i]] <- RunUMAP(sample.list[[i]], dims=1:30)
#  print("umap") 
  #sce.list[[i]] <- SCTransform(sce.list[[i]])
  #sample.list[[i]] <- CreateAssayObject(data=sce@assays@data$logcounts)
  #sce.list[[i]] <- logNormCounts(sce.list[[i]])
}

print("before plots")

#sce.list <- merge.list

#sce.list <- sample.list[[1]]
#for (i in 2:length(sample.list)){
#  sce.list <- merge(sce.list, sample.list[[i]])
#}

#DefaultAssay(sce.list) <- "SCT"
#print("default SCT")

#sce.list <- lapply(sce.list,FUN=SCTransform)
#sce.list <- lapply(sce.list, FUN=RunPCA)
#sce.list <- lapply(sce.list, FUN=RunUMAP)

#print("before d1 d2")

#d1 <- DimPlot(sce.list[[1]],
#        group.by='SampleID')
#d2 <- Dimplot(sce.list[[1]],
#	group.by='SampleID', assay='RNA') 
#d1+d2
#ggpubr::ggarrange(d1,d2, ncol=2)

#d1 <- DimPlot(sample.list[[1]],
#        group.by='SampleID')
#d2 <- DimPlot(sample.list[[1]],
#        group.by='SampleID')

#ggpubr::ggarrange(d1,d2, ncol=2)

#FeaturePlot(sce.list[[1]], features=c('CD3D', 'CD4', 'CD8A', 'CD44', 'CD68', 'ITGAM', 'CD274', 'CD28', 'PDCD1'), order=T)
#DimPlot(sce.list[[1]],
#        group.by='seurat_clusters', label=T)

#sample.list <- sce.list

print("after plots")

for (i in 1:length(sample.list)){
	DefaultAssay(sample.list[[i]]) <- "RNA"
}
print("default RNA")
#saveRDS(sce.list, file=paste(filepath,'outputs/sample_list.RDS',sep=""))

####### INTEGRATION #######

#samples.seurat <- as.Seurat(sample.list)

print(sample.list[[1]])

samples <- sample.list[[1]]
for (i in 2:length(sample.list)){
	samples <- merge(samples, sample.list[[i]])
}

print("finished merge")

adata <- convertFormat(samples, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)

#print(scvi$data)

#adata <- adata.list[[i]]
#for (i in 2:length(adata.list)){
#	adata <- merge(adata, adata.list[[i]])
#}
#adata <- concat(*adata.list)

print(adata)

print("finished convert")

print("starting integration...")
scvi$model$SCVI$setup_anndata(adata, batch_key='orig.ident')
print("setup anndata")
model=scvi$model$SCVI(adata)
model$train()

print("model trained")
saveRDS(adata, file=paste(filepath,'outputs/scvi_adata.RDS',sep=""))
saveRDS(model, file=paste(filepath,'outputs/scvi_model.RDS',sep=""))

latent=model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(samples)
samples[["scvi"]] <- CreateDimReducObject(embeddings=latent, key="scvi_", assay=DefaultAssay(samples))

print("after integration")
saveRDS(samples, file=paste(filepath,'outputs/scvi_samples.RDS',sep=""))

samples <- FindNeighbors(samples, reduction="scvi", dims=1:30)
samples <- FindClusters(samples, resolution=0.2)
samples <- RunUMAP(samples, dims=1:30, reduction="scvi", n.components=2)

merge.immune <- samples

print("before dimplot")

DimPlot(merge.immune,
        group.by='SampleID')
FeaturePlot(merge.immune, features=c('CD3D', 'CD4', 'CD8A', 'CD44', 'CD68', 'ITGAM', 'CD274', 'CD28', 'PDCD1'), order=T)
DimPlot(merge.immune,
        group.by='seurat_clusters', label=T)

dev.off()

saveRDS(merge.immune, file=paste(filepath,'outputs/scvi_immune.RDS',sep=""))
