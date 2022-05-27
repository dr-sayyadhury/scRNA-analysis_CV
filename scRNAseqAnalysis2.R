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

# ```
# 
# ### subsettting out immune cells and cancr cells
# ```{r }

filepath <- '~/projects/def-gbader/cvolk/'

m <- readRDS(file=paste(filepath,'outputs/m.RDS',sep=""))


`%ni%` <- Negate(`%in%`)
immune <- subset(m, subset=seurat_clusters  %in% c(0,1,8,9))
tumor <- subset(m, subset=seurat_clusters  %ni% c(0,1,8,9))
FeaturePlot(immune, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
FeaturePlot(tumor, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
DimPlot(immune, 
        group.by='seurat_clusters', label=T)
DimPlot(tumor, 
        group.by='seurat_clusters', label=T)

# ```
# ### split sample/sce list into tumour & immune cells ### 
# ```{r split, include=F}
sample.list.immune <- SplitObject(immune, split.by = "orig.ident")
sce.list.immune <- lapply(sample.list.immune, as.SingleCellExperiment)

sample.list.tumor <- SplitObject(tumor, split.by = "orig.ident")
sce.list.tumor <- lapply(sample.list.tumor, as.SingleCellExperiment)

# ```
# 
# ### transform again - immune & tumour cells ###
# ```{r transform, include=F}
for (i in 1:length(sample.list.immune)) {
  sample.list.immune[[i]] <- SCTransform(sample.list.immune[[i]])
  #sample.list[[i]] <- NormalizeData(sample.list[[i]]) %>% FindVariableFeatures(., selection.method = 'vst', nfeatures=3000) %>% ScaleData(.)
  #clust <- quickCluster(sce.list[[i]])
  #sce <- computeSumFactors(sce.list[[i]], cluster=clust, min.mean=0.1)
  #sce <- scater::logNormCounts(sce)
  #sample.list[[i]] <- CreateAssayObject(data=sce@assays@data$logcounts)
}
for (i in 1:length(sample.list.tumor)) {
  sample.list.tumor[[i]] <- SCTransform(sample.list.tumor[[i]])
}
# ```
# 
# ### merge tumour & immune ### 
# ```{r merge, include=F }

m.immune <- sample.list.immune[[1]]
for (i in 2:length(sample.list.immune)){
  m.immune <- merge(m, sample.list.immune[[i]])
}

m.tumor <- sample.list.tumor[[1]]
for (i in 2:length(sample.list.tumor)){
  m.tumor <- merge(m, sample.list.tumor[[i]])
}
# ```
# 
# ### transform merged lists ###
# ```{r, include=F }

m.immune <- SCTransform(m.immune)
m.tumor <- SCTransform(m.tumor)

m.immune <- RunPCA (m.immune)
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

m.immune <- FindClusters(m.immune)
m.tumor <- FindClusters(m.tumor)

# ```
# 
# ### Plotting reclustered immune & tumor###
# ```{r plot}

pdf(file=paste(filepath,'outputs/cluster_plots_2.pdf',sep=""), width=15, height=9)

DimPlot(m.immune, 
        group.by='SampleID')
FeaturePlot(m.immune, features=c('PTPRC', 'CD24', 'ITGAM', 'HLA-A'), order=T)
DimPlot(m.immune, 
        group.by='seurat_clusters', label=T)

DimPlot(m.tumor, 
        group.by='SampleID')
FeaturePlot(m.tumor, features=c('PTPRC', 'CD24', 'ITGAM', 'HLA-A'), order=T)
DimPlot(m.tumor, 
        group.by='seurat_clusters', label=T)

dev.off()

saveRDS(m.immune, file=paste(filepath,'outputs/m_immune.RDS',sep=""))
saveRDS(m.tumor, file=paste(filepath,'outputs/m_tumor.RDS',sep=""))

