# Packages and Work Dir

setwd("~/Projects/ImmuOmics_XL")

library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(future)

sample_names <- basename(list.dirs(path = "./Matrix", recursive = FALSE))

# Read matrix
obj_list <- lapply(sample_names, function(s) {
  cat(paste0("Start to read the matrix of ", s, "\n"))
  
  mat <- Read10X(data.dir = file.path("Matrix", s))
  obj <- CreateSeuratObject(counts = mat, project = "ImmuneMM", min.cells = 3)
  obj$sample <- s
  obj$group <- gsub("[0-9]", "", s)
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = grep("^Hba|^Hbb", rownames(obj), value = TRUE))
  
  cat(paste0("Read the matrix of ", s, "\n"))
  
  return(obj)
})

names(obj_list) <- sample_names

for (i in 1:length(obj_list)) {
  cat("===", names(obj_list)[i], "===\n")
  cat("cell number BEFORE filter:", ncol(obj_list[[i]]), "\n")
}
rm(i)

# Merge: confirm criteria
obj_merge <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = "ImmuneMM_merged"
)
obj_merge$sample <- factor(obj_merge$sample, levels = names(obj_list))

VlnPlot(obj_merge, features = "percent.mt", group.by = "sample", pt.size = 0) + NoLegend()
VlnPlot(obj_merge, features = "percent.mt", group.by = "sample") + NoLegend()
VlnPlot(obj_merge, features = "percent.HB", group.by = "sample", pt.size = 0) + NoLegend()
VlnPlot(obj_merge, features = "percent.HB", group.by = "sample") + NoLegend()
VlnPlot(obj_merge, features = "nFeature_RNA", group.by = "sample", pt.size = 0) + NoLegend()
VlnPlot(obj_merge, features = "nFeature_RNA", group.by = "sample") + NoLegend()
VlnPlot(obj_merge, features = "nCount_RNA", group.by = "sample", pt.size = 0) + NoLegend()
VlnPlot(obj_merge, features = "nCount_RNA", group.by = "sample") + NoLegend()


# QC
obj_list <- lapply(obj_list, function(obj) {
  subset(
    obj,
    subset =
      nFeature_RNA > 250 &
      nFeature_RNA < 5500 &
      percent.mt < 25 
  )
})

# QC statistics
for (i in 1:length(obj_list)) {
  cat("===", names(obj_list)[i], "===\n")
  cat("cell number BEFORE filter:", ncol(obj_list[[i]]), "\n")
}
rm(i)

# SCTransform
# obj_list <- lapply(obj_list, function(obj) {
#   SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
# })
for (i in seq_along(obj_list)) {
  s <- names(obj_list)[i]
  cat("===", s, "===\n")
  
  obj_list[[i]] <- SCTransform(
    obj_list[[i]],
    vars.to.regress = "percent.mt",
    conserve.memory = TRUE,     
    verbose = FALSE
  )
  gc()
}
rm(i, s)


obj_list <- lapply(obj_list, function(obj) {
  DefaultAssay(obj) <- "SCT"
  obj
})

# Integration 
dims_int <- 1:30

features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)

# [太大了不好运行]
# anchors <- FindIntegrationAnchors(
#   object.list = obj_list,
#   normalization.method = "SCT",
#   anchor.features = features,
#   dims = dims_int
# )

obj_list <- lapply(obj_list, function(obj) {
  obj <- RunPCA(obj, verbose = FALSE)
  obj
})

# [改成rpca的方法还是太大没法运行]
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features,
  dims = dims_int,
  reduction = "rpca"
)

obj_int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = dims_int)


# Harmony
obj_merge <- obj_merge <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = "ImmuneMM_merged"
)

