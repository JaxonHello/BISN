# ==================================
# Cell Cluster Annotation
# ==================================

# ========== Library ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(RColorBrewer)
})

# ========== Configuration & Initiation ==========
# == work dir == 
config_work_dir <- "~/Projects/ImmuOmics_XL/V01"
setwd(config_work_dir)
# == data: Seurat obj ==
config_obj <- ""
obj <- readRDS(config_obj)
load("/Users/jaxonhe/Projects/ImmuOmics_XL/V01/sc_obj_0108.RData")

rm(config_work_dir, config_obj); gc()

# ========== Function ==========
DimPlot(obj, reduction = "umap", group.by = "MainType")

# subtype
maintype <- "NK-like T cell"
subobj <- subset(obj, subset = MainType == maintype)
p <- DimPlot(subobj, reduction = "umap", group.by = "SubType")
p
ggsave(filename = paste0("~/Desktop/0_annotation/", maintype, ".png"),
       plot = p, width = 6, height = 5, dpi = 300)

rm(subobj, p); gc()

# c("Pro-B", "Pre-B", "Cycling Pre-B", "Mature B cells", "APC+ Naive B")
# c("Promyelocyte", "Proliferating Promyelocyte", "Transitional Myelocyte", 
# "Mature Activated Neutrophil", "Activated Neutrophil")
# c("cDC1", "cDC2", "Cycling Monocyte", "IFN-stimulated myeloid", 
# "Inflammatory Monocyte", "Monocyte-derived Macrophages")
# c("Resident-like Macrophage")
# c("Dendritic Cell")
# c("NK Cell")
# c("CD4+ T Cell", "CD8+ T Cell", "ILC2", "NK-like T Cell", "γδ T Cell")
# c("Plasma Cell")
# c("Erythrocytes")
# c("Mast Cell Precursor")


subtype_levels <- c(
  "Pro-B", "Pre-B", "Cycling Pre-B", "Mature B cells", "APC+ Naive B",
  "Promyelocyte", "Proliferating Promyelocyte", "Transitional Myelocyte",
  "Mature Activated Neutrophil", "Activated Neutrophil",
  "cDC1", "cDC2", "Cycling Monocyte", "IFN-stimulated myeloid",
  "Inflammatory Monocyte", "Monocyte-derived Macrophages",
  "Resident-like Macrophage",
  "Dendritic Cell",
  "NK Cell",
  "CD4+ T Cell", "CD8+ T Cell", "ILC2", "NK-like T Cell", "γδ T Cell",
  "Plasma Cell",
  "Erythrocytes",
  "Mast Cell Precursor"
)

obj$SubType <- factor(as.character(obj$SubType), levels = subtype_levels)
setdiff(unique(as.character(obj$SubType)), subtype_levels)
rm(subtype_levels)


## --- B cell internal ---
ProB_marker <- c("Dntt","Rag1","Rag2")
PreB_marker <- c("Vpreb3","Cd24a")
Cycling_PreB_marker <- c("Mki67","Top2a","Tpx2","Ccnb2")
Mature_B_cells_marker <- c("Ms4a1","Bank1","Cd22")
APC_Naive_B_marker <- c("H2-Aa","H2-Ab1","H2-Eb1","Cd74", "Ccr7")

## --- Neutrophil / granulocyte internal ---
Promyelocyte_marker <- c("Mpo","Elane","Prtn3","Ctsg","Azurocidin1")
Proliferating_Promyelocyte_marker <- c("Mki67","Top2a","Tpx2","Ube2c","Cenpf")
Transitional_Myelocyte_marker <- c("Camp","Ltf","Mmp8")
Activated_Neutrophil_marker <- c("Il1b","Cxcl2","Nfkbia","Lcn2")
Mature_Activated_Neutrophil_marker <- c("Cxcr2","Retnlg")

## --- Monocyte/Macrophage internal ---
Cycling_Monocyte_marker <- c("Mki67","Top2a","Tpx2","Ube2c","Cenpf","Ccnb1")
IFN_stimulated_myeloid_marker <- c("Isg15","Ifit1","Ifit3","Irf7","Stat1","Rsad2","Oas1a","Cxcl10")
Inflammatory_Monocyte_marker <- c("Ly6c2","Ccr2","S100a8","S100a9","Il1b")
Monocyte_derived_Macrophages_marker <- c("Apoe","Lgals3","Spp1","Lpl","Ctsd", "Tyrobp")
Resident_like_Macrophage_marker <- c("Mertk","Timd4","Lyve1","Folr2","Cd163")

## --- DC internal ---
cDC1_marker <- c("Xcr1","Clec9a","Irf8","Wdfy4")
cDC2_marker <- c("Cd209a","Clec10a","Irf4")
Dendritic_Cell_marker <- c("Flt3","Cst3","Ccr7","Ltb")

## --- T/NK/ILC internal ---
CD4_T_marker <- c("Cd4","Il7r","Ltb","Ccr7")
CD8_T_marker <- c("Cd8a","Cd8b1","Gzmk","Cxcr3")
NK_Cell_marker <- c("Nkg7","Klrk1","Klrd1","Prf1","Gzmb","Tyrobp","Fcerg1")
NK_like_T_Cell_marker <- c("Nkg7","Klrk1","Gzmb","Prf1")
gd_T_marker <- c("Trdc", "Trgc2","Rorc","Ccr6")
ILC2_marker <- c("Gata3","Rora","Il1rl1","Areg")

## --- Other single classes ---
Plasma_Cell_marker <- c("Jchain","Xbp1","Mzb1","Sdc1","Prdm1","Derl3")
Erythrocytes_marker <- c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2","Slc4a1","Klf1")
Mast_Cell_Precursor_marker <- c("Kit","Gata2","Cpa3","Tpsb2","Ms4a2","Hdc")

marker_list <- list(
  ProB = ProB_marker,
  PreB = PreB_marker,
  Cycling_PreB = Cycling_PreB_marker,
  Mature_B = Mature_B_cells_marker,
  APC_Naive_B = APC_Naive_B_marker,
  
  Promyelocyte = Promyelocyte_marker,
  Prolif_Promyelocyte = Proliferating_Promyelocyte_marker,
  Transitional_Myelocyte = Transitional_Myelocyte_marker,
  Activated_Neutrophil = Activated_Neutrophil_marker,
  Mature_Activated_Neutrophil = Mature_Activated_Neutrophil_marker,
  
  Cycling_Monocyte = Cycling_Monocyte_marker,
  IFN_myeloid = IFN_stimulated_myeloid_marker,
  Inflammatory_Monocyte = Inflammatory_Monocyte_marker,
  MonoDerived_Macro = Monocyte_derived_Macrophages_marker,
  ResidentLike_Macro = Resident_like_Macrophage_marker,
  
  cDC1 = cDC1_marker,
  cDC2 = cDC2_marker,
  DC = Dendritic_Cell_marker,
  
  CD4_T = CD4_T_marker,
  CD8_T = CD8_T_marker,
  NK = NK_Cell_marker,
  NK_like_T = NK_like_T_Cell_marker,
  gd_T = gd_T_marker,
  ILC2 = ILC2_marker,
  
  Plasma = Plasma_Cell_marker,
  Erythrocytes = Erythrocytes_marker,
  Mast_Precursor = Mast_Cell_Precursor_marker
)

all_markers <- unique(unlist(marker_list, use.names = FALSE))

all_markers <- all_markers[all_markers %in% rownames(obj)]

## ========= DotPlot =========
DotPlot(
  object = obj,
  features = all_markers,
  group.by = "SubType"
) + RotatedAxis()


# B marker
B_markers <- unique(c(
  ProB_marker,
  PreB_marker,
  Cycling_PreB_marker,
  Mature_B_cells_marker,
  APC_Naive_B_marker
))

p <- DotPlot(
  object = subset(obj, subset = MainType == "B cell"),
  features = B_markers,
  group.by = "SubType"
) + RotatedAxis()
ggsave(filename = paste0("~/Desktop/0_annotation/Bcell_sub_marker.png"),
       plot = p, width = 12, height = 7, dpi = 300)

# --- Granulocyte internal ---
Gran_markers <- unique(c(
  Promyelocyte_marker,
  Proliferating_Promyelocyte_marker,
  Transitional_Myelocyte_marker,
  Activated_Neutrophil_marker,
  Mature_Activated_Neutrophil_marker
))

Gran <- subset(obj, subset = MainType == "Granulocyte")
rm(obj); gc()

p <- DotPlot(
  object = subset(obj, subset = MainType == "Granulocyte"),
  features = Gran_markers,
  group.by = "SubType"
) + RotatedAxis()
p
ggsave(filename = paste0("~/Desktop/0_annotation/Gran_sub_marker.png"),
       plot = p, width = 12, height = 7, dpi = 300)

# --- Monocyte/Macrophage internal ---
Mono_markers <- unique(c(
  cDC1_marker,
  cDC2_marker,
  Cycling_Monocyte_marker,
  IFN_stimulated_myeloid_marker,
  Inflammatory_Monocyte_marker,
  Monocyte_derived_Macrophages_marker,
  Resident_like_Macrophage_marker
))

p <- DotPlot(
  object = subset(obj, subset = MainType %in% c("Macrophage", "Monocyte")),
  features = Mono_markers,
  group.by = "SubType"
) + RotatedAxis()
p
ggsave(filename = paste0("~/Desktop/0_annotation/Mono_sub_marker.png"),
       plot = p, width = 12, height = 7, dpi = 300)

# --- NK-like T internal ---
NKlikeT_markers <- unique(c(
  NK_Cell_marker, 
  CD4_T_marker,  
  CD8_T_marker,  
  ILC2 = ILC2_marker,
  NK_like_T_Cell_marker,
  gd_T_marker    
))

p <- DotPlot(
  object = subset(obj, subset = MainType %in% c("NK-like T cell", "NK cell")),
  features = NKlikeT_markers,
  group.by = "SubType"
) + RotatedAxis()
p
ggsave(filename = paste0("~/Desktop/0_annotation/NKT_sub_marker.png"),
       plot = p, width = 12, height = 7, dpi = 300)



