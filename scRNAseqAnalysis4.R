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

m.immune <- readRDS(file=paste(filepath,'outputs/m_immune_final.RDS',sep=""))
m.tumor <- readRDS(file=paste(filepath,'outputs/m_tumor_final.RDS',sep=""))

# ```
# ### RE-SUBSET CELLS ###
# ```{r }
`%ni%` <- Negate(`%in%`)
immune1 <- subset(m.immune, subset=seurat_clusters  %ni% c(19)) # everything except these clusters from prev. clustered immune cells
immune2 <- subset(m.tumor, subset=seurat_clusters %in% c(8)) # these clusters from prev. clustered tumor cells
#immune <- merge(immune1, immune2)

tumor1 <- subset(m.immune, subset=seurat_clusters  %in% c(19)) # these clusters from prev. clustered immune cells
tumor2 <- subset(m.tumor, subset=seurat_clusters %ni% c(8)) # everything except these clusters from prev. clustered tumor cells
#tumor <- merge(tumor1, tumor2)


FeaturePlot(immune1, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
FeaturePlot(immune2, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
FeaturePlot(tumor1, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
FeaturePlot(tumor2, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)
DimPlot(immune1, 
        group.by='seurat_clusters', label=T)
DimPlot(immune2, 
        group.by='seurat_clusters', label=T)
DimPlot(tumor1, 
        group.by='seurat_clusters', label=T)
DimPlot(tumor2, 
        group.by='seurat_clusters', label=T)

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
# ```
# ### cluster ###
# ```{r cluster}
# m.immune <- FindClusters(m.immune, resolution = 0.2)
# m.tumor <- FindClusters(m.tumor, resolution = 0.8)
# 
# ```

### Plotting reclustered immune & tumor###
# ```{r plot}

pdf(file=paste(filepath,'outputs/cluster_plots_4.pdf',sep=""), width=15, height=9)

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

saveRDS(m.immune, file=paste(filepath,'outputs/m_immune_final_2.RDS',sep=""))
saveRDS(m.tumor, file=paste(filepath,'outputs/m_tumor_final_2.RDS',sep=""))
