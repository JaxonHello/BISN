# run_3_integrations_min.R
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

# 1) load obj_list
cat("=== Start Load Data ===\n")
load(RDATA_FILE)  
stopifnot(exists("obj_list"))
cat("=== Load Data Finish ===\n")

# 2) basic setup: ensure SCT default + add sample + unique cell names
if (is.null(names(obj_list)) || any(names(obj_list) == "")) {
  names(obj_list) <- paste0("sample_", seq_along(obj_list))
  cat("the name has changed to sample_index\n")
}
obj_list <- mapply(function(obj, s) {
  DefaultAssay(obj) <- "SCT"
  obj$sample <- s
  obj <- RenameCells(obj, add.cell.id = s)  # 避免不同样本 barcode 重名
  obj
}, obj_list, names(obj_list), SIMPLIFY = FALSE)

# =========================
# C) Harmony (embedding correction)
# =========================
cat("===== C) Harmony (embedding correction) =====\n")
cat("merge obj...\n")
obj_merge <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = paste0(OUT_PREFIX, "_merged")
)
cat("finish merge obj\n")

DefaultAssay(obj_merge) <- "SCT"
obj_merge <- RunPCA(obj_merge, npcs = max(dims_int), verbose = FALSE)

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

# =========================
# A) Seurat default integration (SCT)
# =========================
cat("===== A) Seurat default integration (SCT) =====\n")
features <- SelectIntegrationFeatures(obj_list, nfeatures = nfeatures_int)
obj_list_prep <- PrepSCTIntegration(obj_list, anchor.features = features)

cat("find anchors...\n")
anchors_default <- FindIntegrationAnchors(
  object.list = obj_list_prep,
  normalization.method = "SCT",
  anchor.features = features,
  dims = dims_int
)
cat("finish find anchors\n")

cat("intergrate data...\n")
obj_int_default <- IntegrateData(
  anchorset = anchors_default,
  normalization.method = "SCT",
  dims = dims_int
)
cat("finish intergrate data\n")

cat("cluster and UMAP...\n")
DefaultAssay(obj_int_default) <- "integrated"
obj_int_default <- RunPCA(obj_int_default, npcs = max(dims_int), verbose = FALSE)
obj_int_default <- RunUMAP(obj_int_default, reduction = "pca", dims = dims_int)
obj_int_default <- FindNeighbors(obj_int_default, reduction = "pca", dims = dims_int)
obj_int_default <- FindClusters(obj_int_default, resolution = resolution)
cat("finish cluster and UMAP\n")

cat("save RDS...\n")
saveRDS(obj_int_default, paste0(OUT_PREFIX, "_integrated_seurat_default.rds"))
cat("finish save RDS\n")

rm(features, obj_list_prep, anchors_default, obj_int_default); gc()
cat("==========\n")

# =========================
# B) Seurat RPCA integration (SCT + RPCA)
# =========================
cat("===== B) Seurat RPCA integration (SCT + RPCA) =====\n")
features <- SelectIntegrationFeatures(obj_list, nfeatures = nfeatures_int)
obj_list_prep <- PrepSCTIntegration(obj_list, anchor.features = features)

obj_list_pca <- lapply(obj_list_prep, function(obj) {
  RunPCA(obj, npcs = max(dims_int), verbose = FALSE)
})

cat("run PCA...\n")
anchors_rpca <- FindIntegrationAnchors(
  object.list = obj_list_pca,
  normalization.method = "SCT",
  anchor.features = features,
  dims = dims_int,
  reduction = "rpca"
)
cat("finish PCA\n")

cat("integrate data...\n")
obj_int_rpca <- IntegrateData(
  anchorset = anchors_rpca,
  normalization.method = "SCT",
  dims = dims_int
)
cat("finish integrate\n")

cat("cluster and UMAP...\n")
DefaultAssay(obj_int_rpca) <- "integrated"
obj_int_rpca <- RunPCA(obj_int_rpca, npcs = max(dims_int), verbose = FALSE)
obj_int_rpca <- RunUMAP(obj_int_rpca, reduction = "pca", dims = dims_int)
obj_int_rpca <- FindNeighbors(obj_int_rpca, reduction = "pca", dims = dims_int)
obj_int_rpca <- FindClusters(obj_int_rpca, resolution = resolution)
cat("finish cluster and UMAP\n")
cat("save RDS...\n")
saveRDS(obj_int_rpca, paste0(OUT_PREFIX, "_integrated_seurat_rpca.rds"))
cat("finish save RDS\n")

rm(features, obj_list_prep, obj_list_pca, anchors_rpca, obj_int_rpca); gc()
cat("==========\n")

cat("DONE.\n")
