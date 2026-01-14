# ==================================
# Cell Composition
# ==================================

# ========== Library ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(rlang)
  library(patchwork)
  library(openxlsx)
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
cell_stackplot <- function(meta.data, group = "seurat_clusters", stack = "sample", color = NULL) {
  stopifnot(is.data.frame(meta.data))
  stopifnot(all(c(group, stack) %in% colnames(meta.data)))
  
  proportion_df <- meta.data %>%
    select(all_of(c(group, stack))) %>%
    group_by(across(all_of(c(group, stack)))) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(across(all_of(group))) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  p <- ggplot(proportion_df, aes(x = .data[[group]], y = proportion, fill = .data[[stack]])) +
    geom_col(position = "stack", width = 0.7, col = "black") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    theme_classic() +
    labs(x = group, y = "Proportion (%)", fill = stack) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  if (!is.null(color)) {
    stopifnot(is.character(color))
    if (is.null(names(color))) stop("color must be a *named* vector: names are levels of `stack`.")
    p <- p + scale_fill_manual(values = color, drop = FALSE)
  }
  
  return(p)
}

# ========== Analysis ==========

# MainType in sample
p <- cell_stackplot(obj@meta.data, group = "sample", stack = "MainType", color = pal_main)
ggsave(filename = "V02/CCom/MainType_in_sample.png", plot = p, width = 6, height = 5, dpi = 300)

# SubType in sample
p <- cell_stackplot(obj@meta.data, group = "sample", stack = "SubType")
ggsave(filename = "V02/CCom/SubType_in_sample.png", plot = p, width = 8, height = 5, dpi = 300)

# SubTyep in MainType in sample
p <- cell_stackplot(obj@meta.data[obj@meta.data$MainType == "B cell", ], group = "sample", stack = "SubType")
p
ggsave(filename = "V02/CCom/SubType_in_Bcell_in_sample.png", plot = p, width = 6, height = 5, dpi = 300)

p <- cell_stackplot(obj@meta.data[obj@meta.data$MainType == "Granulocyte", ], group = "sample", stack = "SubType")
p
ggsave(filename = "V02/CCom/SubType_in_Granulocyte_in_sample.png", plot = p, width = 8, height = 5, dpi = 300)

p <- cell_stackplot(obj@meta.data[obj@meta.data$MainType == "NK-like T cell", ], group = "sample", stack = "SubType")
p
ggsave(filename = "V02/CCom/SubType_in_Tcell_in_sample.png", plot = p, width = 6, height = 5, dpi = 300)

p <- cell_stackplot(obj@meta.data[obj@meta.data$MainType == "Monocyte", ], group = "sample", stack = "SubType")
p
ggsave(filename = "V02/CCom/SubType_in_Monocyte_in_sample.png", plot = p, width = 8, height = 5, dpi = 300)

# MainType in group

# SubType in group

# SubType in MainType in group

# CellCycle in sample
p <- cell_stackplot(obj@meta.data, group = "sample", stack = "S.Score")
p
ggsave(filename = "V02/CCom/CellCycle_in_sample.png", plot = p, width = 6, height = 5, dpi = 300)


# CellCycle in MainType in sample
DimPlot(obj, reduction = "umap", group.by = "S.Score")

for (maintype in unique(obj@meta.data$MainType)) {
  print(paste0("========", maintype, "========"))
  p <- cell_stackplot(obj@meta.data[obj@meta.data$MainType == maintype, ], group = "sample", stack = "S.Score")
  p
  ggsave(filename = paste0("V02/CCom/CellCycle_in_", maintype, "_in_sample.png"), plot = p, width = 8, height = 5, dpi = 300)
}

# CellCycl in SubTyep in sample
for (subtype in unique(obj@meta.data$SubType)) {
  print(paste0("========", subtype, "========"))
  p <- cell_stackplot(obj@meta.data[obj@meta.data$SubType == subtype, ], group = "sample", stack = "S.Score")
  p
  ggsave(filename = paste0("V02/CCom/CellCycle_in_", subtype, "_in_sample.png"), plot = p, width = 8, height = 5, dpi = 300)
}



