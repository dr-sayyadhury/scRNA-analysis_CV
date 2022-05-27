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

filepath <- '~/projects/def-gbader/cvolk/'

m.immune <- readRDS(file=paste(filepath,'outputs/m_immune.RDS',sep=""))
m.tumor <- readRDS(file=paste(filepath,'outputs/m_tumor.RDS',sep=""))

# ```
# ### RE-SUBSET CELLS ###
# ```{r }
`%ni%` <- Negate(`%in%`)
immune1 <- subset(m.immune, subset=seurat_clusters  %in% c(0,2,3,4,6,9,12,14,16,17,20,22,23)) # these clusters from prev. clustered immune cells
immune2 <- subset(m.tumor, subset=seurat_clusters %in% c(1,2,3,4,5,9,10,12,18,19,21,22)) # these clusters from prev. clustered tumor cells
#immune <- merge(immune1, immune2)

tumor1 <- subset(m.immune, subset=seurat_clusters  %ni% c(0,2,3,4,6,9,12,14,16,17,20,22,23)) # everything except these clusters from prev. clustered immune cells
tumor2 <- subset(m.tumor, subset=seurat_clusters %ni% c(1,2,3,4,5,9,10,12,18,19,21,22)) # everything except these clusters from prev. clustered tumor cells
#tumor <- merge(tumor1, tumor2)


#FeaturePlot(immune1, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
#FeaturePlot(immune2, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
#FeaturePlot(tumor1, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
#FeaturePlot(tumor2, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
#DimPlot(immune1, 
#        group.by='seurat_clusters', label=T)
#DimPlot(immune2, 
#        group.by='seurat_clusters', label=T)
#DimPlot(tumor1, 
#        group.by='seurat_clusters', label=T)
#DimPlot(tumor2, 
#        group.by='seurat_clusters', label=T)

immune1_counts <- immune1@assays$RNA@counts
immune2_counts <- immune2@assays$RNA@counts

tumor1_counts <- tumor1@assays$RNA@counts
tumor2_counts <- tumor2@assays$RNA@counts

immune1 <- CreateSeuratObject(counts = immune1_counts, project="GBM", meta.data=immune1[[]])
immune2 <- CreateSeuratObject(counts = immune2_counts, project="GBM", meta.data=immune2[[]])

tumor1 <- CreateSeuratObject(counts = tumor1_counts, project="GBM", meta.data=tumor1[[]])
tumor2 <- CreateSeuratObject(counts = tumor2_counts, project="GBM", meta.data=tumor2[[]])


#sample.list.immune1 <- SplitObject(immune1, split.by = "orig.ident")
#sample.list.immune2 <- SplitObject(immune2, split.by = "orig.ident")
#sce.list.immune1 <- lapply(sample.list.immune1, as.SingleCellExperiment)
#sce.list.immune2 <- lapply(sample.list.immune2, as.SingleCellExperiment)

#sample.list.tumor1 <- SplitObject(tumor1, split.by = "orig.ident")
#sample.list.tumor2 <- SplitObject(tumor2, split.by = "orig.ident")
#sce.list.tumor1 <- lapply(sample.list.tumor1, as.SingleCellExperiment)
#sce.list.tumor2 <- lapply(sample.list.tumor2, as.SingleCellExperiment)

#for (i in 1:length(sample.list.immune1)) {
 # sample.list.immune1[[i]] <- SCTransform(sample.list.immune1[[i]])
  #sample.list[[i]] <- NormalizeData(sample.list[[i]]) %>% FindVariableFeatures(., selection.method = 'vst', nfeatures=3000) %>% ScaleData(.)
  #clust <- quickCluster(sce.list[[i]])
  #sce <- computeSumFactors(sce.list[[i]], cluster=clust, min.mean=0.1)
  #sce <- scater::logNormCounts(sce)
  #sample.list[[i]] <- CreateAssayObject(data=sce@assays@data$logcounts)
#}

#for (i in 1:length(sample.list.immune2)) {
 # sample.list.immune2[[i]] <- SCTransform(sample.list.immune2[[i]])
#}

#for (i in 1:length(sample.list.tumor1)) {
 # sample.list.tumor1[[i]] <- SCTransform(sample.list.tumor1[[i]])
#}

#for (i in 1:length(sample.list.tumor2)) {
 # sample.list.tumor2[[i]] <- SCTransform(sample.list.tumor2[[i]])
#}

#m.immune <- immune1[[1]]
#for (i in 2:length(immune1)){
#  m.immune <- merge(m.immune, sample.list.immune1[[i]])
#}
#for (i in 2:length(sample.list.immune2)){
#  m.immune <- merge(m.immune, sample.list.immune2[[i]])
#}

#m.tumor <- sample.list.tumor1[[1]]
#for (i in 2:length(sample.list.tumor1)){
#  m.tumor <- merge(m.tumor, sample.list.tumor1[[i]])
#}
#for (i in 2:length(sample.list.tumor2)){
#  m.tumor <- merge(m.tumor, sample.list.tumor2[[i]])
#}

# ```
# ### MERGE ###
# ```{r merge. include=F}
m.immune <- merge(immune1, immune2)
m.tumor <- merge(tumor1, tumor2)
# ```
# 
# ### transform merged lists ###
# ```{r, include=F }

m.immune <- SCTransform(m.immune)
m.tumor <- SCTransform(m.tumor)

m.immune <- RunPCA(m.immune)
m.tumor <- RunPCA(m.tumor)

m.immune <- RunUMAP(m.immune, 
                    dims=1:30,
                    n.neighbors = 10)
m.tumor <- RunUMAP(m.tumor, 
                   dims=1:30,
                   n.neighbors = 10)

m.immune <- FindNeighbors(m.immune, 
                          reduction ="pca")
m.tumor <- FindNeighbors(m.tumor, 
                         reduction ="pca")
# ```
# ### cluster ###
# ```{r cluster}
m.immune <- FindClusters(m.immune, resolution = 0.8)
m.tumor <- FindClusters(m.tumor, resolution = 0.8)
# 
# ```

### Plotting reclustered immune & tumor###
# ```{r plot}

pdf(file=paste(filepath,'outputs/cluster_plots_3.pdf',sep=""), width=15, height=9)

DimPlot(m.immune, 
        group.by='SampleID')
FeaturePlot(m.immune, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A', 'PECAM1'), order=T)
DimPlot(m.immune, 
        group.by='seurat_clusters', label=T)

DimPlot(m.tumor, 
        group.by='SampleID')
FeaturePlot(m.tumor, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A', 'PECAM1'), order=T)
DimPlot(m.tumor, 
        group.by='seurat_clusters', label=T)

dev.off()

saveRDS(m.immune, file=paste(filepath,'outputs/m_immune_final.RDS',sep=""))
saveRDS(m.tumor, file=paste(filepath,'outputs/m_tumor_final.RDS',sep=""))
