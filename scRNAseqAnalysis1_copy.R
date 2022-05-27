# ---
#   title: "Untitled"
# output: html_document
# date: '2022-05-06'
# ---
#   
#   ```{r setup, include=FALSE}
# knitr::opts_chunk$set(echo = TRUE)
# ```


# Packages
# ```{r packages, include=FALSE}

library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)
library(gridExtra)
library(ggpubr)
library(scran)


filepath <- '~/projects/def-gbader/cvolk/'

m <- readRDS(file=paste(filepath,'outputs/m.RDS',sep=""))

###--------------- DATA LOADING -------------------- ###
# ```{r read.data, include=FALSE}
#counts <- read.csv(file=paste(filepath,'single_cell_dataset/SCP503/other/Richards_NatureCancer_GBM_scRNAseq_counts.csv',sep=""),
#                   header = TRUE,
#                   row.names = 1,
#                   sep = ",")

#colnames(counts) <- gsub(x = colnames(counts),
#                         pattern = "\\.",
#                         replacement = "-")

#meta <- read.csv(file=paste(filepath,'single_cell_dataset/SCP503/other/Richards_NatureCancer_GBM_scRNAseq_meta.csv',sep=""),
#                 header = TRUE,
#                 row.names = 1,
#                 sep = ",")

#all(rownames(meta)==colnames(counts))

# ```

### SEURAT PIPELINE ###
# ```{r seurat.obj, include=TRUE}
#m <- CreateSeuratObject(counts = counts,
#                          project = "GBM",
#                          meta.data = meta)


#m <- subset(GBM, subset=SampleID %in% c("G1003-A_T", "G910-A_T","G945-I_T", "G946-I_T","G967-A_T", "G983-A_T"))

# ```

### PLOTTING ###
# ```{r plots.QC, fig.height=4, fig.width=9}

#m <- SCTransform(m)
#m <- RunPCA (m)
#m <- RunUMAP(m, 
#             dims=1:30,
#             n.neighbors = 10)

#m <- FindNeighbors(m, 
#                   reduction ="pca")


#for(i in seq(0,1,0.1)){
#       m <- FindClusters(m, resolution = i)
#}


#pdf(file=paste(filepath,'outputs/clustree_1.pdf', sep=""), width=15, height=9)
#clustree(m, prefix="SCT_snn_res.") +
#       ggtitle("Clustering tree (m)")

#dev.off()

m  <- FindClusters(m, resolution=0.2)

# ```
# ### PLOT CLUSTERS ###
# ```{r }

pdf(file=paste(filepath,'outputs/plots_1.pdf',sep=""), width=15, height=9)
DimPlot(m, 
        group.by='SampleID')

FeaturePlot(m, features=c('PTPRC', 'CD68', 'SOX2', 'HLA-A'), order=T)

DimPlot(m, 
        group.by='seurat_clusters', label=T)
dev.off()

saveRDS(m, file=paste(filepath,'outputs/m.RDS',sep=""))

### STOP EXECUTION HERE - SAVE IN RDS ###
