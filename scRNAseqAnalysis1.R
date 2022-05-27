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
#library(colorspace)
#library(STACAS)
library(clustree)
library(gridExtra)
library(ggpubr)
#library(egg)
#library(scales)
library(scran)
#library(effectsize)

#```

filepath <- '~/projects/def-gbader/cvolk/'

###--------------- DATA LOADING -------------------- ###
# ```{r read.data, include=FALSE}
counts <- read.csv(file=paste(filepath,'single_cell_dataset/SCP503/other/Richards_NatureCancer_GBM_scRNAseq_counts.csv',sep=""),
                   header = TRUE,
                   row.names = 1,
                   sep = ",")

colnames(counts) <- gsub(x = colnames(counts),
                         pattern = "\\.",
                         replacement = "-")

meta <- read.csv(file=paste(filepath,'single_cell_dataset/SCP503/other/Richards_NatureCancer_GBM_scRNAseq_meta.csv',sep=""),
                 header = TRUE,
                 row.names = 1,
                 sep = ",")

all(rownames(meta)==colnames(counts))

# ```

### SEURAT PIPELINE ###
# ```{r seurat.obj, include=TRUE}
GBM <- CreateSeuratObject(counts = counts,
                          project = "GBM",
                          meta.data = meta)

#GBM[["ribo"]] <- PercentageFeatureSet(GBM, pattern = "^RP[SL]")
#GBM[["mito"]] <- PercentageFeatureSet(GBM, pattern = "^MT-")

#saveRDS(GBM, paste(filepath,'outputs/GBM.rds',sep="")
# ```
### SAMPLE LISTS ###
# ```{r }
#red_list <- list()
#GBM <- readRDS(GBM, paste(filepath,'outputs/GBM.rds',sep="")
#red_list <- c('')

#GBMr <- subset(GBM, subset=orig.ident %in% red_list)
#GBMr <- GBM

#rm(GBM)
#sample.list <- SplitObject(GBMr, split.by = "orig.ident")
#sce.list <- lapply(sample.list, as.SingleCellExperiment)

GBMr <- subset(GBM, subset=SampleID %in% c(“G1003-A_T”, “G910-A_T”,“G945-I_T”, “G946-I_T”,“G967-A_T”, “G983-A_T”))

# ```

### PLOTTING ###
# ```{r plots.QC, fig.height=4, fig.width=9}
ggplot(GBMr[[]], aes(nCount_RNA, mito, color=ribo)) +
  geom_jitter(size=0.3) +
  theme_pubr() +
  #scale_y_log10() +
  scale_x_log10()
#  scale_color_continuous_sequential(palette='ag_Sunset', rev=T) 
# facet_wrap(~orig.ident, nrow=2)

ggplot(GBMr[[]], aes(nGene, mito, color=PatientID)) +
  geom_jitter(size=0.3) +
  theme_pubr() +
  #scale_y_log10() +
  scale_x_log10()
  #scale_color_continuous_sequential(palette='ag_Sunset', rev=T) +
  facet_wrap(~orig.ident, nrow=2)

ggplot(GBMr[[]], aes(nCount_RNA, nGene, color=ribo)) +
  geom_jitter(size=0.3) +
  theme_pubr() +
  scale_y_log10() +
  scale_x_log10()
 # scale_color_continuous_sequential(palette='ag_Sunset', rev=T)

ggplot(GBMr[[]], aes(nCount_RNA, nGene, color=mito)) +
  geom_jitter(size=0.3) +
  theme_pubr() +
  scale_y_log10() +
  scale_x_log10()
 # scale_color_continuous_sequential(palette='ag_Sunset', rev=T)

ggplot(GBMr[[]], aes(nCount_RNA, nGene, color=orig.ident)) +
  geom_jitter(size=0.3) +
  theme_pubr() +
  scale_y_log10() +
  scale_x_log10() 

#ggplot(GBMr[[]], aes(nCount_RNA, nGene, color=nGene<200)) +
 # geom_jitter(size=0.3) +
 # theme_pubr() +
 # scale_y_log10() +
 # scale_x_log10() 

# ```

### INITIAL INTEGRATION ###
# ``` {r integration, include=TRUE}

#for (i in 1:length(sample.list)) {
 # sample.list[[i]] <- SCTransform(sample.list[[i]])
  #sample.list[[i]] <- NormalizeData(sample.list[[i]]) %>% FindVariableFeatures(., selection.method = 'vst', nfeatures=3000) %>% ScaleData(.)
  #clust <- quickCluster(sce.list[[i]])
  #sce <- computeSumFactors(sce.list[[i]], cluster=clust, min.mean=0.1)
  #sce <- scater::logNormCounts(sce)
  #sample.list[[i]] <- CreateAssayObject(data=sce@assays@data$logcounts)
#}
# ```
# 
# ### MERGING ###
# ```{r, include=F }

#m <- sample.list[[1]]
#for (i in 2:length(sample.list)){
#  m <- merge(m, sample.list[[i]])
#}
# ```
# 
# ### CLUSTERING ###
# ```{r, include=F }

m <- SCTransform(m)
m <- RunPCA (m)
m <- RunUMAP(m, 
             dims=1:30,
             n.neighbors = 10)

m <- FindNeighbors(m, 
                   reduction ="pca")

m  <- FindClusters(m)

# ```
# ### PLOT CLUSTERS ###
# ```{r }

pdf(file=paste(filepath,'outputs/cluster_plots_1.pdf',sep=""), width=15, height=9)
DimPlot(m, 
        group.by='SampleID')

FeaturePlot(m, features=c('PTPRC', 'CD24', 'FAM107A', 'HLA-A'), order=T)

DimPlot(m, 
        group.by='seurat_clusters', label=T)
dev.off()

saveRDS(m, file=paste(filepath,'outputs/m.RDS',sep=""))

### STOP EXECUTION HERE - SAVE IN RDS ###
