# Packages and Work Dir

setwd("~/Projects/ImmuOmics_XL")

library(Seurat)
library(DoubletFinder)
library(dplyr)

sample_name = "AD1"

# Read Data
data <- Read10X(data.dir = paste0("Matrix/", sample_name))
obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3)
rm(data)

obj$sample <- sample_name
obj$group <- sample_name

# QC
# 1) percent.mt: 线粒体基因占比
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
# 2) percent.HB: 血红蛋白基因占比，用来评估红细胞污染
obj[["percent.HB"]] <- PercentageFeatureSet(obj, features =  grep("^Hba|^Hbb", rownames(obj), value = TRUE))
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")) # pt.size = 0
print(paste0("Cell Number before QC:", ncol(obj)))

obj <- subset(obj, subset = nFeature_RNA > 250 & nFeature_RNA < 7500 & percent.mt < 25)
# obj <- subset(obj, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 25 & percent.mt < 10)
print(paste0("Cell Number after QC:", ncol(obj)))

##### Basicinfo_nGene_nUMI_mito.txt & QC criteria #####
if(FALSE) {
  # Compare with company result
  ref <- read.table(file = "Basicinfo_nGene_nUMI_mito.txt", header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
  ref <- ref[grepl(paste0("^", sample_name), ref$Cell), ]
  print(paste0("Cell Number before QC by 10x company:", nrow(ref)))
  
  meta.table <- obj@meta.data
  meta.table$Cell <- paste0(sample_name, "_", rownames(meta.table))
  
  addi <- meta.table[!meta.table$Cell %in% ref$Cell, ] # find the reason
  # ! ref remove all the cells with percent.mt > 80
  
  ref <- read.csv("AD1.qc_filter_annotation_table.csv")
  # max percent.mt = 25, don't care about the percent.HB
  
  rm(addi, meta.table, ref)
}
#####

obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, npcs = 100)
ElbowPlot(obj, ndims = 60)

dim_use <- 50
obj <- FindNeighbors(obj, dims = 1:dim_use)
obj <- FindClusters(obj, resolution = 0.8)
obj <- RunUMAP(obj, dims = 1:dim_use)

DimPlot(obj, label = TRUE, repel = TRUE, reduction = "umap")

##### Doublet Find [No need] ####
if(FALSE) {
  set.seed(2512)
  
  sweep.res  <- paramSweep(obj, PCs = 1:dim_use, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  bcmvn2 <- subset(bcmvn, pK >= 0.01)
  best_pK <- bcmvn2$pK[which.max(bcmvn2$BCmetric)]
  best_pK
  
  pN = 0.25
  nExp <- round(ncol(obj) * 0.055)
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp_adj <- round(nExp_use * (1 - homotypic.prop))
  
  obj <- doubletFinder_v3(
    obj,
    PCs = dims_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_adj,
    reuse.pANN = FALSE,
    sct = TRUE
  )
  
  meta.data <- obj@meta.data
  pann_col <- grep("^pANN_", colnames(meta.data), value = TRUE)
  dfcol_col <- grep("^DF\\.classifications_", colnames(meta.data), value = TRUE)
  
  
}
#####

