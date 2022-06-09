
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

filepath <- '~/projects/def-gbader/cvolk/'

m.immune <- readRDS(file=paste(filepath,'outputs/m_immune.RDS',sep=""))

### Initial cluster visualization ###

pdf(file=paste(filepath,'outputs/immune_plots.pdf',sep=""), width=15, height=9)

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

### Re-normalization & Integration ###

for (i in 1:length(sample.list)) {
  sample.list[[i]] <- SCTransform(sample.list[[i]])
}

integ.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures=3000)
m.immune <- PrepSCTIntegration(object.list = sample.list, anchor.features = integ.features)
integ.anchors <- FindIntegrationAnchors(object.list = m.immune, normalization.method = "SCT", anchor.features = integ.features)
integ.immune <- IntegrateData(anchorset = integ.anchors, normalization.method = "SCT", k.weight=35)

### Visualization of integrated data ###

#merge.immune <- integ.immune[[1]]
#for (i in 2:length(integ.immune)){
#  merge.immune <- merge(merge.immune, integ.immune[[i]])
#}
merge.immune <- integ.immune

merge.immune <- RunPCA(merge.immune)
merge.immune <- RunUMAP(merge.immune, dims=1:30, n.neighbours=10, reduction="pca")

DimPlot(merge.immune,
        group.by='SampleID')

dev.off()

saveRDS(merge.immune, file=paste(filepath,'outputs/integ_immune.RDS',sep=""))
