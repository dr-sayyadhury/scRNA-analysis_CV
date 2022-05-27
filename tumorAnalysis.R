
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

m.tumor <- readRDS(file=paste(filepath,'outputs/m_tumor.RDS',sep=""))



pdf(file=paste(filepath,'outputs/tumor_plots.pdf',sep=""), width=15, height=9)

DimPlot(m.tumor,
        group.by='SampleID')
FeaturePlot(m.tumor, features=c('CD3D', 'CD4', 'CD8A', 'CD44', 'CD68', 'ITGAM', 'CD274', 'CD28', 'PDCD1'), order=T)
DimPlot(m.tumor,
        group.by='seurat_clusters', label=T)

#DimPlot(m.tumor,
#        group.by='SampleID')
#FeaturePlot(m.tumor, features=c('PTPRC', 'CD24', 'ITGAM', 'HLA-A'), order=T)
#DimPlot(m.tumor,
#        group.by='seurat_clusters', label=T)

dev.off()

