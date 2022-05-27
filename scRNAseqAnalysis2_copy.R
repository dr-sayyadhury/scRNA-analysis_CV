library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)
library(gridExtra)
library(ggpubr)
library(scran)

# ```
# 
# ### subsettting out immune cells and cancr cells
# ```{r }

filepath <- '~/projects/def-gbader/cvolk/'

m <- readRDS(file=paste(filepath,'outputs/m.RDS',sep=""))


`%ni%` <- Negate(`%in%`)
m.immune <- subset(m, subset=seurat_clusters  %in% c(0,1,4,8))
m.tumor <- subset(m, subset=seurat_clusters  %ni% c(0,1,4,8))

# 
# ### merge tumour & immune ### 
# ```{r merge, include=F }

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

#for(i in seq(0,1,0.1)){
#	m.immune <- FindClusters(m.immune, resolution = i)
#}

#for(i in seq(0,1,0.1)){
#        m.tumor <- FindClusters(m.tumor, resolution = i)
#}

#pdf(file=paste(filepath,'outputs/clustree_2.pdf', sep=""), width=15, height=9)
#clustree(m.immune, prefix="SCT_snn_res.") +
#	ggtitle("Clustering tree (m.immune)")
#clustree(m.tumor, prefix="SCT_snn_res.") +
#        ggtitle("Clustering tree (m.tumor)")

#dev.off()


m.immune <- FindClusters(m.immune, resolution=0.1)
m.tumor <- FindClusters(m.tumor, resolution=0.1)

# ```
# 
# ### Plotting reclustered immune & tumor###
# ```{r plot}

pdf(file=paste(filepath,'outputs/plots_2.pdf',sep=""), width=15, height=9)

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

