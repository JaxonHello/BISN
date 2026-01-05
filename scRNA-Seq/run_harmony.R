# run_harmony.R
setwd("~/ImmuOmics")

library(Seurat)
library(harmony)

# --------- config ----------
RDATA_FILE <- "SCT_obj.RData"
OUT_PREFIX <- "ImmuneMM"
dims_int <- 1:30
nfeatures_int <- 3000
resolution <- 0.8
# --------------------------

cat("=== Start Load Data ===\n")
load(RDATA_FILE)
stopifnot(exists("obj_list"))
cat("=== Load Data Finish ===\n")

if (is.null(names(obj_list)) || any(names(obj_list) == "")) {
  names(obj_list) <- paste0("sample_", seq_along(obj_list))
  cat("the name has changed to sample_index\n")
}

obj_list <- mapply(function(obj, s) {
  DefaultAssay(obj) <- "SCT"
  obj$sample <- s
  obj <- RenameCells(obj, add.cell.id = s)
  obj
}, obj_list, names(obj_list), SIMPLIFY = FALSE)


cat("select integration features...\n")
features_use <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures_int)
cat("features:", length(features_use), "\n")

cat("===== C) Harmony (embedding correction) =====\n")
cat("merge obj...\n")
obj_merge <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  project = paste0(OUT_PREFIX, "_merged")
)
cat("finish merge obj\n")

DefaultAssay(obj_merge) <- "SCT"
VariableFeatures(obj_merge) <- features_use
obj_merge <- RunPCA(
  obj_merge,
  npcs = max(dims_int),
  features = features_use,
  verbose = FALSE
)

cat("run harmony...\n")
obj_merge <- RunHarmony(
  object = obj_merge,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = dims_int,
  assay.use = "SCT"
)
cat("finish run harmony\n")

cat("cluster and UMAP...\n")
obj_merge <- RunUMAP(obj_merge, reduction = "harmony", dims = dims_int)
obj_merge <- FindNeighbors(obj_merge, reduction = "harmony", dims = dims_int)
obj_merge <- FindClusters(obj_merge, resolution = resolution)
cat("finish cluster and UMAP\n")

cat("save RDS...\n")
saveRDS(obj_merge, paste0(OUT_PREFIX, "_harmony.rds"))
cat("finish save RDS\n")
cat("==========\n")