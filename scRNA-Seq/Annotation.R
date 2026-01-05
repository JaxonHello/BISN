# Packages and Work Dir

setwd("~/Projects/ImmuOmics_XL")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(openxlsx)
library(SingleCellExperiment)
library(slingshot)
library(harmony)
library(RColorBrewer)
library(tidyr)
library(CellChat)

# config
dims_int <- 1:30
resolution <- 0.8

# Marker list
# B cell marker
B_marker <- c("Cd79a", "Cd79b", "Cd19", "Ms4a1", "Cd74", "Cd22", "Cd37", "Pax5", "Ebf1", "Pou2af1")
# Erythrocytes(红细胞)，在骨髓中属于是红系细胞，
ET_marker <- c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Gata1", "Klf1", "Epor")
# Hb-是血红细胞的alpha和beta链，Gata1, Klf1, Epor属于红系分化基因

# Granulocytes(粒细胞)：包含中性粒，嗜酸性粒细胞，嗜碱性粒细胞。
GL_marker <- c("Ly6g", "S100a8", "S100a9", "Csf3r", "Mpo", "Elane", "Gfi1", "Cebpe", "Lcn2", "Camp", "Ltf", "Fcnb", "Cd177")

MA_marker <- c("Adgre1", "Csf1r", "Lyz2", "Lgals3", "Fcgr1", "Cd68", "Clec4f", "Vsig4", "Timd4")
MO_marker <- c("Ly6c2", "Ccr2", "Lst1", "S100a8", "S100a9")
NK_marker <- c("Ncr1", "Klrk1", "Klrb1c", "Nkg7", "Prf1", "Gzmb", "Ifng")

# Load combined data
obj <- readRDS("All_sample_combined.rds") # company result
# obj <- readRDS("ImmuneMM_harmony.rds") # our result

# Evaluate the integration performance
DimPlot(obj, reduction = "harmony", group.by = "sample") + ggtitle("Harmony reduction")
DimPlot(obj, reduction = "umap", group.by = "sample") + ggtitle("Harmony UMAP")

obj <- RunUMAP(obj, reduction = "pca", dims = dims_int, reduction.name = "umap_orig")
# obj <- FindNeighbors(obj, reduction = "pca", dims = dims_int, graph.name = c("pcaTest_nn", "pcaTest_snn"))
# obj <- FindClusters(obj, resolution = resolution, cluster.name = "cluster_noInte", graph = "pcaTest_snn")

DimPlot(obj, reduction = "umap_orig", group.by = "sample") + ggtitle("PCA reduction") +
  xlab("UMAP_1") +
  ylab("UMAP_2")

# Find markers
obj <- PrepSCTFindMarkers(object = obj, assay = "SCT")
markers <- FindAllMarkers(obj)

All_markers <- read.csv("allmarkers_annotation.csv")

# Dimplot
DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE)

p_list <- lapply(levels(obj$seurat_clusters), function(cl) {
  DimPlot(subset(obj, seurat_clusters == cl), reduction = "umap") +
    ggtitle(paste0("Cluster ", cl)) +
    coord_cartesian(xlim = c(-12,17), ylim = c(-16,10))
})

p_all <- wrap_plots(p_list, ncol = 6)

ggsave(filename = "base/AllCluster.pdf", p_all, width = 24, height = 20, units = "in")


Idents(obj) <- "seurat_clusters"
cl_levels <- levels(Idents(obj))

emb <- as.data.frame(Embeddings(obj, "umap"))
colnames(emb) <- c("UMAP_1","UMAP_2")
emb$cluster <- factor(Idents(obj), levels = cl_levels)
emb$sample <- obj$sample

pal <- hue_pal()(length(cl_levels)); names(pal) <- cl_levels

p_list <- lapply(cl_levels, function(cl) {
  df <- emb[emb$cluster == cl, , drop = FALSE]
  ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(size = 0.2, color = pal[cl]) +
    coord_cartesian(xlim = c(-12,17), ylim = c(-16,10)) +
    ggtitle(paste0("Cluster ", cl), ) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))
})

p_all <- wrap_plots(p_list, ncol = 6)

ggsave(filename = "base/AllCluster_color.pdf", p_all, width = 24, height = 20, units = "in")

rm(emb, p_list, p_all)
gc()


sample_levels <- unique(emb$sample)
sample_pal <- hue_pal()(length(sample_levels))
names(sample_pal) <- sample_levels

p_list <- lapply(cl_levels, function(cl) {
  df <- emb[emb$cluster == cl, , drop = FALSE]
  ggplot(df, aes(UMAP_1, UMAP_2, color = sample)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = sample_pal) +
    coord_cartesian(xlim = c(-12,17), ylim = c(-16,10)) +
    ggtitle(paste0("Cluster ", cl)) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))
})

p_all <- wrap_plots(p_list, ncol = 6)

ggsave(filename = "base/AllCluster_sampleColor.pdf", p_all, width = 24, height = 20, units = "in")

rm(p_list, p_all)
gc()


# 堆叠柱状图，每个cluster中细胞的数量
cluster_col <- "seurat_clusters"  
sample_col  <- "sample"       
df <- obj@meta.data %>%
  transmute(
    cluster = factor(.data[[cluster_col]]),
    sample  = factor(.data[[sample_col]])
  )

tab <- df %>%
  count(cluster, sample, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_count <- ggplot(tab, aes(x = cluster, y = n, fill = sample)) +
  geom_col(width = 0.85) +
  theme_classic() +
  labs(x = "Cluster", y = "Number of cells", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

p_prop <- ggplot(tab, aes(x = cluster, y = prop, fill = sample)) +
  geom_col(width = 0.85, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "Cluster", y = "Composition within cluster", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# p_count
okabe_ito <- c(
  "#E69F00","#56B4E9","#009E73","#F0E442",
  "#0072B2","#D55E00","#CC79A7","#000000"
)

p_prop +
  scale_fill_manual(values = okabe_ito)

# Marker
# 1. 红细胞系：cluster 18似乎只有一部分是，还需要继续做聚类。
FeaturePlot(obj, features = "percent.HB", reduction = "umap")
VlnPlot(obj, features = "percent.HB")
FeaturePlot(obj, features = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"), reduction = "umap")
p <- VlnPlot(obj, features = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"),  ncol = 4)
ggsave(filename = "base/Annotation/Hb_vlnplot.pdf", p, width = 20, height = 5, units = "in")
FeaturePlot(obj, features = c("Gata1", "Klf1", "Epor"), reduction = "umap") # 红系分化基因高表达
p <- VlnPlot(obj, features = c("Gata1", "Klf1", "Epor"),  ncol = 4)
ggsave(filename = "base/Annotation/Hb_vlnplot_v2.pdf", p, width = 20, height = 5, units = "in")

# 2. cluster 17: Ig Gene高表达
# Ig gene expression: 基本上确定为B细胞。
ig_genes <- grep("^Igh|^Igk|^Igl", rownames(obj), value = TRUE)
obj$percent_Ig <- PercentageFeatureSet(obj, pattern = "^Igh|^Igk|^Igl")
FeaturePlot(obj, features = "percent_Ig") 
VlnPlot(obj, features = "percent_Ig")

# 同时看他的marker基因，还有Jchain，这是plasma cell的重要marker，
# 编码的 J 链参与多聚体免疫球蛋白（IgA、IgM）的组装和分泌，几乎只在浆细胞 / 黏膜 B 细胞中高表达。
FeaturePlot(obj, features = "Jchain")
# 但是还有一些别的信号，暂时标注为plasma cell

# 3. B cell
FeaturePlot(obj, features = B_marker)
# Cd79a/b是BCR复合物的组成部分，所有B细胞都会表达，但是在Naive中Cd79a相对来说高于Cd79b
# Ms4a1是成熟B细胞的marker，可以多关注这个，说明6里面应该有一些事未成熟B细胞，Scd1也是一样的功能。
FeaturePlot(obj, features = c("Cd79a", "Cd79b", "Ms4a1", "Cd22", "Bank1", "Fcer2a"))
# Igd基本上可以确定是Naive B cell
# H2-类基因，是小鼠的MHC II类分子，是抗原呈递的核心分子，高表达说明细胞能高效向CD4+T细胞呈递抗原，高表达于cluster8
FeaturePlot(obj, features = c("Ighd", "Cd74", "H2-Aa", "H2-Eb1", "H2-Ab1", "H2-Ob"))
# 除了Scd1之外，基本上都在8群中表达。
FeaturePlot(obj, features = c("Scd1", "Fcmr", "Bank1", "Fcer2a", "Ccr7", "Cxcr5"))
# Bank1说明BCR信号通路处于活跃状态，这和抗原呈递高表达的基因高表达正好吻合。
# Fcmr是IgM Fc受体，Fcer2a是IgM/IgD的辅助受体，基本上确定是成熟且是Naive B细胞
# Ccr7表达于Naive B，Cxcr5说明该亚群参与生发中心的形成，向淋巴滤泡迁移。看来确实是Naive B
FeaturePlot(obj, features = c("Cd74")) # 高水平（辅助 MHC II 类分子装载抗原）

# 3.1 cluster 8: 这部分细胞是成熟初始 B 细胞（Naive B Cell），且具有专职抗原呈递表型，高抗原呈递的Naive B细胞亚群。

# pro/pre-B("Rag1", "Dntt", "Vpreb1", "Vpreb3", "Igll1", "Ebf1") # 这里的Vpreb3不是跟pre-BCR相关
FeaturePlot(obj, features = c("Rag1", "Rag2", "Dntt", "Vpreb1", "Igll1", "Ebf1"))
# 实锤cluster 13是Pro/pre-B了
# Ebf1的特征是这样的，在pre-B和pro-B中高表达，在成熟Naive B细胞中低表达，有一个趋势在。
# Igll1, Dntt, Vpreb1是Pre-B cell receptor，Lef1，Rag1/Rag2也是pre-B的marker
# Vpreb1 + Igll1 + pre-BCR 是surrogate light chain的组成成分，Dntt发生在VDJ重排

FeaturePlot(obj, features = c("Igkc", "Il2ra")) # 这两个一般在pre-B的后期表达，
# pre-B（尤其 late pre-B）**更常见：Il2ra (CD25)、Rag1/2（轻链重排期）、有时 Vpreb1/Igll1 仍可见

# 跟细胞周期，细胞增殖相关的基因【说明下面有一些细胞也正在细胞增殖】
FeaturePlot(obj, features = c("Mki67", "Top2a", "Stmn1"))
FeaturePlot(obj, features = c("Cd79a", "Pax5", "Ebf1", "Rag1", "Igkc", "Il2ra"))
# 基本上实锤了10群是cycling pre-B/或者说large pre-B吧
# 那么13就是pro-B了

FeaturePlot(obj, features = c("Sox4", "Myb", "Ly6d", "Polm", "Fcrla", "Ms4a1"))
FeaturePlot(obj, features = c("Cd93", "Tnfrsf13c"))
# cluster 6总体来说就是B细胞，但是左边是later pre-B，中间是Immature B, 最右边是成熟B细胞，
FeaturePlot(obj, features = c("Bank1", "Fcmr", "Ms4a1", "Ighd")) # 快成熟细胞的maker

# 16似乎跟旁边的表达都很相似，只不过表达都很淡


FeaturePlot(obj, features = c("Igkc", "Iglc3")) # 正好说明小鼠的IgK基因是dominant的
FeaturePlot(obj, features = c("Cd3d", "Cd3e", "Trac", "Nkg7", "Ccl5", "Xcl1", "Cd8a")) 
# cluster 11
# 这群细胞是 T 细胞（而不是 NK），并且偏 细胞毒性/效应型（cytotoxic/effector）T，
# 很可能是 CD8 T（或 NK-like CD8 T）

# cluster 24
FeaturePlot(obj, features = "nCount_RNA")
VlnPlot(obj, features = "nFeature_RNA")

Bcell <- subset(obj, idents = c(13, 10, 16, 24, 6, 25, 8))
saveRDS(Bcell, file = "base/Bcell.rds")
rm(Bcell); gc()

Bmarker_25 <- FindMarkers(obj, 
                          ident.1 = 25, 
                          ident.2 = c(13, 10, 16, 24, 6, 8),
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.5)


# 粒细胞(大概在UMAP的左边，并且大量)
FeaturePlot(obj, features = GL_marker, ncol = 6)


FeaturePlot(obj, features = c("Mpo", "Elane", "Prtn3", "Ctsg", "Mmp8")) # Mmp8不保守
# Camp, Ngp: 4,3 ; Mpo, Elane, Prtn3, Cts: 12, 7的一部分。
# 4，3的分裂活动比12, 7更高。
FeaturePlot(obj, features = c("Mki67", "Top2a", "Stmn1"))
VlnPlot(obj, features = c("Mki67", "Top2a", "Stmn1"))

FeaturePlot(obj, features = c("Ltf", "Camp", "Ngp", "S100a8","S100a9", "Ly6g"))
VlnPlot(obj, features = c("Ltf", "Camp", "Ngp", "Ly6g"), pt.size = 0)

# 高表达Mmp9的中性粒细胞，更易从骨髓迁出到外周炎症部位
# 更偏向于三级
FeaturePlot(obj, features = c("Mmp9"))
VlnPlot(obj, features = c("Mmp9"), pt.size = 0)

# 趋化因子分泌基因，炎症微环境招募，动员免疫细胞，因此离免疫细胞更近
# 可以做调控因子网络分析，细胞通讯。
# Il1b和Ptgs2是中性粒细胞被激活的标记。
FeaturePlot(obj, features = c("Cxcl2", "Ccl6", "Ccrl2", "Il1b", "Ptgs2"))
VlnPlot(obj, features = c("Cxcl2", "Ccl6", "Ccrl2", "Il1b", "Ptgs2"), pt.size = 0)

# 跟病原体应答和中性粒细胞发挥功能，避免凋亡，提供能量等有关
# 说明中性粒细胞正在发挥作用
FeaturePlot(obj, features = c("Nlrp3", "Clec4d", "Clec7a", "Acod1", "Hcar2", "Smox"))
# 这两个并不能作为monocyte的marker，因为在粒细胞中也有表达。
# 并且Ly6c2正好在活化的中性粒细胞中低表达
FeaturePlot(obj, features = c("Lyz2", "Ly6c2"))

#
FeaturePlot(obj, features = c("Ccr2", "Fn1", "F13a1", "Msr1", "Ms4a4c", "Ly86"))

FeaturePlot(obj, features = c("Ms4a7", "FcεRIa", "Tpsab1", "Prss34"))
FeaturePlot(obj, features = c("Ccr2", "Ms4a7", "Ctss", "Ccr2", "Cd19"))

# 以下说的是Macrophage的marker，除了Ms4a7(这个基因比较特异性高表达)
# 单核细胞（Monocyte）和巨噬细胞（Macrophage）在转录组上通常很像，
# 尤其在骨髓/炎症条件下，它们经常处在连续谱系（monocyte → macrophage）
# 区分思路("单核更偏迁移/炎症前体；巨噬更偏吞噬/组织驻留/补体脂代谢")
FeaturePlot(obj, features = c("Csf1r", "Mafb", "Cx3cr1", "Apoe", "Lpl", "Cd68"))
# 以下是巨噬细胞核心标志
FeaturePlot(obj, features = c("Csf1r", "C1qa", "Mafb", "Cd68"))
# 脂代谢最相关基因 and 组织驻留细胞
FeaturePlot(obj, features = c("Apoe", "Lpl", "Cx3cr1"))
# 特异性高表达，肯跟跟亚群有关
FeaturePlot(obj, features = c("Axl", "Ms4a7"))
# 最后注释为组织驻留巨噬细胞 resident-like Macrophage
# 最后确认，排除monocyte污染
FeaturePlot(obj, features = c("Ly6c2", "Ccr2", "S100a8", "S100a9"))

# 经典Monocytes: Ccr2+Ly6c2+Csf1r-C1qa-Mafb-
FeaturePlot(obj, features = c("Ccr2", "Ly6c2", "Csf1r", "C1qa", "Mafb"))

FeaturePlot(subset(obj, !is.na(pANN_0.25_0.09_2438)), 
            features = "pANN_0.25_0.09_2438", cols = c("lightgray", "red"))
VlnPlot(subset(obj, !is.na(pANN_0.25_0.28_2120)), 
            features = "pANN_0.25_0.28_2120", pt.size = 0)

# 先把7-20-5-23注释为Intermediate mono-mac

#
VlnPlot(obj, features = "percent.mt", pt.size = 0)

obj$highlight <- ifelse(obj$seurat_clusters == 24, 
                        paste("Cluster", 24), 
                        "Other Clusters")

FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "highlight",
               cols = c("#C075A6", "lightgrey"), pt.size = 3) + ggtitle("cluster 24")


cl_levels <- levels(Idents(obj))
pal <- setNames(Seurat:::DiscretePalette(length(cl_levels), palette = "polychrome"),
                cl_levels)

df0 <- obj@meta.data[, c("nCount_RNA", "percent.mt"), drop = FALSE]
df0$cluster <- as.character(Idents(obj))

set.seed(1)
bg_max <- 30000  
df_bg <- if (nrow(df0) > bg_max) df0[sample.int(nrow(df0), bg_max), ] else df0

make_scatter_one <- function(cl) {
  df_hl <- df0[df0$cluster == cl, , drop = FALSE]  
  df_bg$grp <- "Other"
  df_hl$grp <- cl
  
  df_plot <- rbind(df_bg, df_hl)
  
  ggplot(df_plot, aes(x = nCount_RNA, y = percent.mt, color = grp)) +
    geom_point(size = 0.25, alpha = 0.6) +
    scale_color_manual(values = c("Other" = "lightgray", cl = pal[as.integer(cl)])) +
    labs(title = paste0("Cluster ", cl)) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

p_list <- lapply(cl_levels, make_scatter_one)

p_all <- wrap_plots(p_list, ncol = 6)

ggsave("base/QC_scatter_by_cluster.png", p_all, width = 24, height = 18, dpi = 300)


# cluster 17
FeaturePlot(obj, features = c("Irtf8", "Flt3", "Tcf4", "Siglech", "Bst2"))
# Mpeg1 在Mon-Mac中也是高表达，说明里面有细胞正在分化为DC
# cluster 17为DC，更准确地说，应该是Plasmacytoid Dendritic Cell, pDC

# cluster 19
# 前两个是肥大细胞的谱系marker，后面跟细胞增殖有关，代表了一种活跃的增殖状态。
FeaturePlot(obj, features = c("Cpa3", "Gata2", "Cdk6", "Ccnd2"))

FeaturePlot(obj, features = c("Ms4a1", "Ighd", "Ighm"))

FeaturePlot(subset(obj, subset = seurat_clusters == 16), features = "percent.mt")
FeaturePlot(subset(obj, subset = seurat_clusters == 16), features = "nCount_RNA")

# 看起来分布不是很奇怪，再细分。
FeaturePlot(subset(obj, subset = seurat_clusters == 16 & percent.mt < 15 & nCount_RNA > 3000), 
            features = "nCount_RNA")


rm(df_bg, df0, p_list, bg_max, cl_levels, GL_marker, i, p_all, pal, make_scatter_one)
gc()
# mainType annotation
anno_table <- read.xlsx("Cell_Annotation.xlsx") %>% select(seurat_clusters, MainType)
anno_table$seurat_clusters <- as.character(anno_table$seurat_clusters)
obj$MainType <- anno_table$MainType[ match(as.character(obj$seurat_clusters), anno_table$seurat_clusters) ]
table(is.na(obj$MainType))
DimPlot(obj, reduction = "umap", group.by = "MainType")
saveRDS(obj, "base/IntewithMainType.rds")

# 开始做细分&减法
# 1. remove low quality
obj <- subset(obj, subset = MainType != "Low quality")

# 2. NK细胞的细分
nkcd8 <- subset(obj, MainType == "NT-like CD8T")

DefaultAssay(nkcd8) <- "RNA"
nkcd8 <- SCTransform(nkcd8, vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(nkcd8) <- "SCT"

nkcd8 <- RunPCA(nkcd8, npcs = 30)
nkcd8 <- FindNeighbors(nkcd8, dims = 1:20)
nkcd8 <- FindClusters(nkcd8, resolution = 0.4)
nkcd8 <- RunUMAP(nkcd8, dims = 1:20)
DimPlot(nkcd8, label = TRUE)


DefaultAssay(nkcd8) <- "RNA"

nkcd8 <- NormalizeData(nkcd8)
nkcd8 <- FindVariableFeatures(nkcd8, nfeatures = 2000)
nkcd8 <- ScaleData(nkcd8, vars.to.regress = "percent.mt")
nkcd8 <- RunPCA(nkcd8)

nkcd8 <- FindNeighbors(nkcd8, dims = 1:30)
nkcd8 <- FindClusters(nkcd8, resolution = 0.9)  
nkcd8 <- RunUMAP(nkcd8, dims = 1:30)
DimPlot(nkcd8, label = TRUE, group.by = "sample") # 比较均匀
DimPlot(nkcd8, label = TRUE)


FeaturePlot(nkcd8, features = c("Cd3d", "Cd3e", "Trac", "Nkg7", "Ccl5", "Xcl1", "Cd8a", "Cd4",
                                "Cd8b1"))
# FeaturePlot(nkcd8, features = c("Ccr7", "Cd3e", "Cd28", "Cd62L"))

# Highlight
emb0 <- as.data.frame(Embeddings(obj, "umap"))
cells_pick <- rownames(emb0)[ Idents(obj)[rownames(emb0)] == "11" & emb0$UMAP_1 < 2 ]

DimPlot(nkcd8, reduction = "umap", cells.highlight = cells_pick) # cluster 2
marker2 <- FindMarkers(nkcd8, ident.1 = 2)
# 2并不是单纯的中性粒细胞
marker5 <- FindMarkers(nkcd8, ident.1 = 5)
# γδ T Cell作为5群细胞
# CD4+ T Cell作为6群
# 7和1是NK细胞
# 8, 3, 4是CD8+ T Cell
# 0, 9, 10, 11, 13, 2是NK-like T Cell
# 12 ILC2

emb0 <- as.data.frame(Embeddings(nkcd8, "umap"))
cells_pick <- rownames(emb0)[ Idents(nkcd8)[rownames(emb0)] == "2"]
DimPlot(obj, reduction = "umap", cells.highlight = cells_pick)

marker0 <- FindMarkers(nkcd8, ident.1 = 0)
marker_12 <- FindMarkers(
  object = nkcd8,
  ident.1 = 12,   
  min.pct = 0.1,
  only.pos = TRUE
)
# 12 ILC2

# 开始注释
emb0 <- as.data.frame(Embeddings(obj, "umap"))
emb0$cell <- rownames(emb0)
cells_gran <- emb0$cell[ Idents(obj)[emb0$cell] == "11" & emb0$UMAP_1 < 2 ]

nkcd8$cellType <- NA_character_

# 按你给的规则赋值（注意：cluster 编号请用字符）
nkcd8$cellType[nkcd8$seurat_clusters %in% c("5")]  <- "γδ T Cell"
nkcd8$cellType[nkcd8$seurat_clusters %in% c("6")]  <- "CD4+ T Cell"
nkcd8$cellType[nkcd8$seurat_clusters %in% c("7","1")] <- "NK Cell"
nkcd8$cellType[nkcd8$seurat_clusters %in% c("8","3","4")] <- "CD8+ T Cell"
nkcd8$cellType[nkcd8$seurat_clusters %in% c("0","9","10","11","13","2", "14")] <- "NK-like T Cell"
nkcd8$cellType[nkcd8$seurat_clusters %in% c("12")] <- "ILC2"

cells_gran_in_nk <- intersect(cells_gran, colnames(nkcd8))
cat("[INFO] cells_gran_in_nkcd8:", length(cells_gran_in_nk), "\n")
nkcd8$cellType[cells_gran_in_nk] <- "Granulocytes"

DimPlot(nkcd8, reduction = "umap", group.by = "cellType")
saveRDS(nkcd8, "Base/nkcd8.rds")

# 
obj$SubType <- NA_character_
obj$SubType[colnames(nkcd8)] <- nkcd8$cellType

DimPlot(
  obj,
  reduction = "umap",
  group.by = "SubType"
)
rm(emb0, marker_12, marker0, marker5, nkcd8, cells_gran, cells_gran_in_nk, cells_pick); gc()


# 将游离点归入粒细胞
emb <- as.data.frame(Embeddings(obj, "umap"))
emb$cell <- rownames(emb)
emb$cluster <- as.character(Idents(obj)[emb$cell])

cells_gran_rule <- emb$cell[
  (emb$cluster == "18" & emb$UMAP_1 < 1) |
    (emb$cluster == "16" & emb$UMAP_1 < 4) |
    (emb$cluster == "9"  & emb$UMAP_1 < 1)
]

cat("[INFO] selected cells:", length(cells_gran_rule), "\n")
obj$MainType[cells_gran_rule] <- "Granulocytes"

DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  cells.highlight = cells_gran_rule,
  cols.highlight = "red",
  sizes.highlight = 1.0
) + ggtitle("Highlight: candidate Granulocytes (by cluster & UMAP_1 threshold)")
# ✅

#########################
# 初步分析
DimPlot(obj, reduction = "umap", group.by = "MainType")
obj$MainType[obj$MainType == "NT-like CD8T"] <- "NK-like T cell"
obj$MainType[obj$SubType == "NK Cell"] <- "NK cell"
DimPlot(obj, reduction = "umap", group.by = "MainType")

obj@meta.data$sample_test <- factor(obj@meta.data$sample, levels = c("AD1", "AD3", "AS1", "AS2", "Aged1", "Control2"))

tab <- obj@meta.data %>%
  transmute(
    sample   = .data[["sample_test"]],
    MainType = as.character(MainType)
  ) %>%
  filter(!is.na(sample), !is.na(MainType)) %>%
  count(sample, MainType, name = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(tab, aes(x = sample, y = prop, fill = MainType)) +
  geom_col(width = 0.85, color = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "Sample", y = "Proportion within sample", fill = "MainType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# 发现AD有更多的resident-like Macrophage，AS1/2几乎不含有NK-like的T cell
# 在Control组中有最多的B cells，自然衰老的次之

# Cell cycle分析
colnames(obj@meta.data)[colnames(obj@meta.data) == "Phase"] <- "Cell_Cycle_Phase"  
colnames(obj@meta.data)[colnames(obj@meta.data) == "S.Score"] <- "S_Phase_Score"   
colnames(obj@meta.data)[colnames(obj@meta.data) == "G2M.Score"] <- "G2M_Phase_Score"  


DimPlot(obj, reduction = "umap", group.by = "Cell_Cycle_Phase", cols = c("G1" = "#3498db", 
                                                              "S" = "#f39c12", 
                                                              "G2M" = "#e74c3c")) +
  ggtitle("Cell Cycle")

# G1期：恢复期，间隙期，S：合成期，G2M: 真正的细胞分裂期

# Cell trajectory分析

# B cell
DimPlot(obj, reduction = "umap", group.by = "MainType")
Bcell <- subset(obj, MainType == "B cell")

DefaultAssay(Bcell) <- "RNA"
Bcell <- NormalizeData(Bcell)
Bcell <- FindVariableFeatures(Bcell, nfeatures = 2000)
Bcell <- ScaleData(Bcell, vars.to.regress = "percent.mt", verbose = FALSE)
Bcell <- RunPCA(Bcell)
Bcell <- FindNeighbors(Bcell, dims = 1:30)
Bcell <- FindClusters(Bcell, resolution = 0.4)
Bcell <- RunUMAP(Bcell, dims = 1:30)

DimPlot(Bcell, reduction = "umap", group.by = "sample")
DimPlot(Bcell, reduction = "umap")


proportion_df <- Bcell@meta.data %>% 
  select(seurat_clusters, sample) %>% 
  dplyr::group_by(seurat_clusters, sample) %>% 
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::mutate(proportion = count / sum(count) * 100)

ggplot(proportion_df, aes(x = seurat_clusters, y = proportion, fill = sample)) +
  geom_col(position = "stack", width = 0.7, col = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "seurat_cluster", y = "Proportion within sample", fill = "MainType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
rm(proportion_df, count_df)

Bcell <- RunHarmony(
  Bcell,
  group.by.vars = "sample",
  dims.use = 1:30
)
Bcell <- RunUMAP(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindNeighbors(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindClusters(Bcell, resolution = 0.8)

DimPlot(Bcell, reduction = "umap", label = TRUE)
DimPlot(Bcell, reduction = "umap", group.by = "sample") # better

# Annotation
DimPlot(Bcell, reduction = "umap", label = TRUE)
DimPlot(Bcell, reduction = "umap")
FeaturePlot(Bcell, features = c("Ighd", "Cd74", "H2-Aa", "H2-Eb1", "H2-Ab1", "H2-Ob"))
VlnPlot(Bcell, features = c("Ighd", "Cd74", "H2-Aa", "H2-Eb1", "H2-Ab1", "H2-Ob"))
# 0, 14, 5的一部分是APC+ Naive B
FeaturePlot(Bcell, features = c("Rag1", "Rag2", "Dntt", "Vpreb1", "Igll1", "Ebf1"))
# 7 10 13(maybe)是Pro-B
FeaturePlot(Bcell, features = c("Igkc", "Il2ra"))
FeaturePlot(Bcell, features = c("Mki67", "Top2a", "Stmn1", "Cd79a", "Pax5", "Ebf1", "Rag1", "Igkc", "Il2ra"))
FeaturePlot(Bcell, features = c("Rag1", "Rag2", "Dntt", "Vpreb1", "Igll1", "Ebf1"))
# 8, 9, 6是Cycling Pre-B
cells <- WhichCells(obj, idents = "24")
cells <- intersect(cells, colnames(Bcell))
DimPlot(Bcell, cells.highlight = cells, reduction = "umap")
# 怀疑这个11群应该是doublets，建议去除，先注释为Doublet.
cells <- WhichCells(Bcell, idents = "11")
DimPlot(Bcell, cells.highlight = cells, reduction = "umap")
# 只有2是Mature B cells，还有5的一部分
FeaturePlot(Bcell, features = c("Ms4a1", "Cd27", "Cd22", "Cd79a", "Il7r"))
# 
marker <- FindMarkers(Bcell, ident.1 = 5, logfc.threshold = 0.5,
                      min.pct = 0.2,
                      only.pos = TRUE)
# 5表达巨噬细胞marker，不像是典型的B细胞，怀疑Doublet，注释Doublet
# 1, 3, 4, 16, 15都注释为pre-B

# 记得加Cd22, Spn = Cd43
FeaturePlot(Bcell, features = c("Ms4a1", "Cd24a", "Spn", "Cd79a", "Il7r", "Ighm"))
VlnPlot(Bcell, features = c("Ms4a1", "Cd24a", "Spn", "Cd79a", "Il7r", "Ighm"), pt.size = 0)

marker <- FindMarkers(Bcell, ident.1 = 12, logfc.threshold = 0.5,
                      min.pct = 0.2,
                      only.pos = TRUE)

# 开始注释：
Bcell$cellType <- "Unannotated"

# 步骤2：按聚类编号逐一赋值（严格遵循你的注释逻辑）
# 0, 14, 5的一部分是APC+ Naive B（注：你后续提到5群表达巨噬细胞marker注释为Doublet，这里优先按你最终结论，若需保留5群部分注释可看后续补充）
Bcell$cellType[Bcell$seurat_clusters %in% c("0", "14")]  <- "APC+ Naive B"
# 7 10 13(maybe)是Pro-B
Bcell$cellType[Bcell$seurat_clusters %in% c("7", "10", "13")]  <- "Pro-B"
# 8, 9, 6是Cycling Pre-B
Bcell$cellType[Bcell$seurat_clusters %in% c("8", "9", "6", "12")]  <- "Cycling Pre-B"
# 11群是Doublets
Bcell$cellType[Bcell$seurat_clusters %in% c("11")]  <- "Doublet"
# 2是Mature B cells
Bcell$cellType[Bcell$seurat_clusters %in% c("2")]  <- "Mature B cells"
# 5群表达巨噬细胞marker，注释为Doublet
Bcell$cellType[Bcell$seurat_clusters %in% c("5")]  <- "Doublet"
# 1, 3, 4, 16, 15都注释为pre-B
Bcell$cellType[Bcell$seurat_clusters %in% c("1", "3", "4", "15", "16")]  <- "Pre-B"

DimPlot(Bcell, reduction = "umap", group.by = "cellType")
saveRDS(Bcell, "Base/Bcell.rds")

# 
Bcell <- subset(Bcell, cellType != "Doublet")
DefaultAssay(Bcell) <- "RNA"
Bcell <- NormalizeData(Bcell)
Bcell <- FindVariableFeatures(Bcell, nfeatures = 2000)
Bcell <- ScaleData(Bcell, vars.to.regress = "percent.mt", verbose = FALSE)
Bcell <- RunPCA(Bcell)

Bcell <- RunHarmony(Bcell, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:30)
Bcell <- RunUMAP(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindNeighbors(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindClusters(Bcell, resolution = 0.8)

DimPlot(Bcell, reduction = "umap", group.by = "cellType") + ggtitle("B cell subtype")

# DimPlot(Bcell, reduction = "pca", group.by = "cellType")
rm(marker, pt, sce, cells, cells_c13, start_cl)

Bcell <- readRDS("Base/Bcell.rds")
obj$SubType[intersect(colnames(obj), colnames(Bcell))] <- Bcell$cellType[intersect(colnames(obj), colnames(Bcell))]
DimPlot(obj, group.by = "SubType")
rm(Bcell) ; gc()


# B细胞轨迹分析
sce <- as.SingleCellExperiment(Bcell)
reducedDim(sce, "HARMONY") <- Embeddings(Bcell, "harmony")
reducedDim(sce, "UMAP") <- Embeddings(Bcell, "umap")

# 找到起点
cells <- WhichCells(Bcell, idents = "16")
cells <- intersect(cells, colnames(obj))
DimPlot(obj, cells.highlight = cells, reduction = "umap")
start_cl <- "7"

# 费时
sce <- slingshot(
  sce,
  clusterLabels = "cellType",
  reducedDim    = "UMAP",
  start.clus    = "Pro-B"
)

pt <- slingPseudotime(sce)

plot(reducedDims(sce)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("Cell trajectory"),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) 
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)] 

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1, main = "B cell trajectory") 
lines(SlingshotDataSet(sce), lwd=2, col="black")
legend("right",
       legend = paste0("lineage",1:6),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)


# filter obj
DimPlot(obj, reduction = "umap", group.by = "MainType")
cells_Bobj <- rownames(obj@meta.data)[obj@meta.data$MainType %in% c("B cell")]

obj$pseudotime <- NA_real_
cells_pt <- intersect(rownames(pt), colnames(obj))
obj$pseudotime[cells_pt] <- as.numeric(pt[cells_pt, "Lineage1"])

rm(pt, cells_pt)
gc()
saveRDS(sce, "Base/Bcell_trajectory.rds")
rm(Bcell, sce); gc()

FeaturePlot(obj, features = "pseudotime")

drop_cells <- rownames(obj@meta.data)[
  obj@meta.data$MainType == "B cell" & is.na(obj@meta.data$pseudotime)
]
keep_cells <- setdiff(colnames(obj), drop_cells)
obj <- subset(obj, cells = keep_cells)
rm(keep_cells, drop_cells, plotcol, colors); gc()


DimPlot(obj, reduction = "umap", group.by = "MainType")

proportion_df <- obj@meta.data %>% 
  filter(MainType == "B cell") %>%
  select(sample, SubType) %>% 
  dplyr::group_by(sample, SubType) %>% 
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(proportion = count / sum(count) * 100)

proportion_df$sample <- factor(proportion_df$sample, levels = c("AD1", "AD3", "AS1", "AS2", "Aged1", "Control2"))

ggplot(proportion_df, aes(x = sample, y = proportion, fill = SubType)) +
  geom_col(position = "stack", width = 0.7, col = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "Sample", y = "Proportion within SubType", fill = "SubType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

rm(proportion_df, count_df)

df <- obj@meta.data %>%
  transmute(
    cell   = rownames(obj@meta.data),
    sample = as.character(.data[["sample"]]),
    MainType = as.character(MainType),
    pt     = as.numeric(pseudotime)
  ) %>%
  filter(MainType == "B cell", !is.na(pt), !is.na(sample))

binwidth <- 0.5
df_bin <- df %>%
  mutate(
    bin = floor(pt / binwidth) * binwidth,
    bin_mid = bin + binwidth/2
  ) %>%
  dplyr::count(sample, bin_mid, name = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>% 
  ungroup()

ggplot(df_bin, aes(x = bin_mid, y = prop, color = sample)) +
  # geom_col(aes(fill = sample), position = "identity",
  #          alpha = 0.18, width = binwidth, linewidth = 0) +
  # geom_line(linewidth = 0.75) +
  geom_smooth(se = FALSE, method = "loess", span = 0.2, linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(
    x = "B-cell pseudotime (Lineage1, binned)",
    y = "Fraction within each sample (B cells only)",
    color = "Sample",
    fill  = "Sample"
  )

rm(df, df_bin, binwidth)

### Gran ! 跟和B cell同样的方法分析##
DimPlot(obj, group.by = "MainType")
Gran <- subset(obj, MainType == "Granulocytes")

# save memory 
Gran@graphs <- list()
obj@graphs <- list()

DefaultAssay(Gran) <- "RNA"
Gran <- NormalizeData(Gran)
Gran <- FindVariableFeatures(Gran, nfeatures = 2000)
Gran <- ScaleData(Gran, vars.to.regress = "percent.mt", verbose = FALSE)
Gran <- RunPCA(Gran)
Gran <- FindNeighbors(Gran, dims = 1:30)
Gran <- FindClusters(Gran, resolution = 0.4)
Gran <- RunUMAP(Gran, dims = 1:30)

DimPlot(Gran, reduction = "umap", group.by = "sample")

Gran <- RunHarmony(
  Gran,
  group.by.vars = "sample",
  dims.use = 1:30
)
Gran <- RunUMAP(Gran, reduction = "harmony", dims = 1:30)
Gran <- FindNeighbors(Gran, reduction = "harmony", dims = 1:30)
Gran <- FindClusters(Gran, resolution = 0.3)

DimPlot(Gran, reduction = "umap", label = TRUE)
DimPlot(Gran, reduction = "umap", group.by = "sample") # better

# !! Very Mixed cluster

# Annotation

FeaturePlot(Gran, features = c("Mpo", "Elane", "Prtn3", "Ctsg", "Mmp8"), reduction = "umap") # Mmp8不保守
# Camp, Ngp: 4,3 ; Mpo, Elane, Prtn3, Cts: 12, 7的一部分。
# 4，3的分裂活动比12, 7更高。
FeaturePlot(Gran, features = c("Mki67", "Top2a", "Stmn1"))
# VlnPlot(obj, features = c("Mki67", "Top2a", "Stmn1"))

FeaturePlot(Gran, features = c("Ltf", "Camp", "Ngp", "S100a8","S100a9", "Ly6g"))
VlnPlot(Gran, features = c("Ltf", "Camp", "Ngp", "Ly6g"), pt.size = 0)

# 高表达Mmp9的中性粒细胞，更易从骨髓迁出到外周炎症部位
# 更偏向于三级(Mmp9更倾向于晚期)
FeaturePlot(Gran, features = c("Mmp9"))
VlnPlot(Gran, features = c("Mmp9"), pt.size = 0)

# 趋化因子分泌基因，炎症微环境招募，动员免疫细胞，因此离免疫细胞更近
# 可以做调控因子网络分析，细胞通讯。
# Il1b和Ptgs2是中性粒细胞被激活的标记。(Il1b和Ptgs2意味着被激活)
FeaturePlot(Gran, features = c("Cxcl2", "Ccl6", "Ccrl2", "Il1b", "Ptgs2"))
VlnPlot(Gran, features = c("Cxcl2", "Ccl6", "Ccrl2", "Il1b", "Ptgs2"), pt.size = 0)

# 跟病原体应答和中性粒细胞发挥功能，避免凋亡，提供能量等有关
# 说明中性粒细胞正在发挥作用
FeaturePlot(Gran, features = c("Nlrp3", "Clec4d", "Clec7a", "Acod1", "Hcar2", "Smox"))
# 这两个并不能作为monocyte的marker，因为在粒细胞中也有表达。
# 并且Ly6c2正好在活化的中性粒细胞中低表达
FeaturePlot(Gran, features = c("Lyz2", "Ly6c2"))

#
FeaturePlot(Gran, features = c("Ccr2", "Fn1", "F13a1", "Msr1", "Ms4a4c", "Ly86"))
FeaturePlot(Gran, features = c("Ms4a7", "FcεRIa", "Tpsab1", "Prss34"))
FeaturePlot(Gran, features = c("Ccr2", "Ms4a7", "Ctss", "Ccr2", "Cd19"))

# Unmature Gran (promyelocyte) c8
FeaturePlot(Gran, features = c("Mpo", "Elane"))

# Later phase: cluster1 注释为成熟粒细胞
FeaturePlot(Gran, features = c("Mmp9", "Mmp8", "Ly6g", "S100a8", "Ltf"))
VlnPlot(Gran, features = c("Mmp9", "Mmp8", "Ly6g", "S100a8", "Ltf"), pt.size = 0)
# Ly6G⁺Ltf⁺Mmp9⁻ ：Transitional Myelocyte cluster2


# Mature Activated Gran (Activated Marker)
FeaturePlot(Gran, features = c("Il1b", "Ptgs2", "Nlrp3", "Clec4d"))
VlnPlot(Gran, features = c("Il1b", "Ptgs2", "Nlrp3", "Clec4d"), pt.size = 0)
# 0, 7, 10, 3, 9, 6

# Proliferating Promyelocyte cluster 4,5
FeaturePlot(Gran, features = c("Mki67", "Top2a", "Stmn1"))

# Annotation part:
Gran$cellType <- "Unannotated"
Gran$cellType[Gran$seurat_clusters %in% c("8")]  <- "Promyelocyte"
Gran$cellType[Gran$seurat_clusters %in% c("4", "5")]  <- "Proliferating Promyelocyte"
Gran$cellType[Gran$seurat_clusters %in% c("2")]  <- "Transitional Myelocyte"
Gran$cellType[Gran$seurat_clusters %in% c("1")]  <- "Mature Activated Neutrophil"
Gran$cellType[Gran$seurat_clusters %in% c("0", "7", "10", "3", "9", "6")]  <- "Activated Neutrophil"

DimPlot(Gran, reduction = "umap", group.by = "cellType")


# cell trajectory
# Gran轨迹分析
sce <- as.SingleCellExperiment(Gran)
reducedDim(sce, "HARMONY") <- Embeddings(Gran, "harmony")
reducedDim(sce, "UMAP") <- Embeddings(Gran, "umap")

# 费时
sce <- slingshot(
  sce,
  clusterLabels = "cellType",
  reducedDim    = "UMAP",
  start.clus    = "Promyelocyte"
)

pt <- slingPseudotime(sce)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) 
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)] 

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1, main = "Granulocytes trajectory") 
lines(SlingshotDataSet(sce), lwd=2, col="black")


# 将注释结果添加到obj.meta的SubType里面
obj.meta.beat <- obj.meta
obj.meta[colnames(Gran), "SubType"] <- Gran@meta.data[colnames(Gran), "cellType"]
obj.meta[colnames(Gran), "pseudotime"] <- pt

proportion_df <- obj.meta %>% 
  filter(MainType == "Granulocytes") %>%
  select(sample, SubType) %>% 
  dplyr::group_by(sample, SubType) %>% 
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(proportion = count / sum(count) * 100)

proportion_df$sample <- factor(proportion_df$sample, levels = c("AD1", "AD3", "AS1", "AS2", "Aged1", "Control2"))

ggplot(proportion_df, aes(x = sample, y = proportion, fill = SubType)) +
  geom_col(position = "stack", width = 0.7, col = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "Sample", y = "Proportion within SubType", fill = "SubType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


df <- obj.meta %>%
  transmute(
    cell   = rownames(obj.meta),
    sample = as.character(.data[["sample"]]),
    MainType = as.character(MainType),
    pt     = as.numeric(pseudotime)
  ) %>%
  filter(MainType == "Granulocytes", !is.na(pt), !is.na(sample))

binwidth <- 0.5
df_bin <- df %>%
  mutate(
    bin = floor(pt / binwidth) * binwidth,
    bin_mid = bin + binwidth/2
  ) %>%
  dplyr::count(sample, bin_mid, name = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>% 
  ungroup()

ggplot(df_bin, aes(x = bin_mid, y = prop, color = sample)) +
  # geom_col(aes(fill = sample), position = "identity",
  #          alpha = 0.18, width = binwidth, linewidth = 0) +
  # geom_line(linewidth = 0.75) +
  geom_smooth(se = FALSE, method = "loess", span = 0.2, linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(
    x = "B-cell pseudotime (Lineage1, binned)",
    y = "Fraction within each sample (Granulocytes only)",
    color = "Sample",
    fill  = "Sample"
  )
saveRDS(Gran, "Base/Gran.rds")
saveRDS(sce, "Base/Gran_trajectory.rds")

# 最后一堆，单核细胞&巨噬细胞混合物。
DimPlot(obj, reduction = "umap", group.by = "MainType")

MM <- subset(obj, MainType == "Intermediate mono-macro")

DefaultAssay(MM) <- "RNA"
MM <- NormalizeData(MM)
MM <- FindVariableFeatures(MM, nfeatures = 2000)
MM <- ScaleData(MM, vars.to.regress = "percent.mt", verbose = FALSE)
MM <- RunPCA(MM)
MM <- FindNeighbors(MM, dims = 1:30)
MM <- FindClusters(MM, resolution = 0.4)
MM <- RunUMAP(MM, dims = 1:30)

DimPlot(MM, reduction = "umap", group.by = "sample")

MM <- RunHarmony(
  MM,
  group.by.vars = "sample",
  dims.use = 1:30
)
MM <- RunUMAP(MM, reduction = "harmony", dims = 1:30)
MM <- FindNeighbors(MM, reduction = "harmony", dims = 1:30)
MM <- FindClusters(MM, resolution = 0.3)

DimPlot(MM, reduction = "umap", label = TRUE)
DimPlot(MM, reduction = "umap", group.by = "sample")

proportion_df <- MM@meta.data %>% 
  select(seurat_clusters, sample) %>% 
  dplyr::group_by(seurat_clusters, sample) %>% 
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::mutate(proportion = count / sum(count) * 100)

ggplot(proportion_df, aes(x = seurat_clusters, y = proportion, fill = sample)) +
  geom_col(position = "stack", width = 0.7, col = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "Sample", y = "Proportion within SubType", fill = "SubType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))



DimPlot(MM, reduction = "umap", label = TRUE)
DimPlot(MM, reduction = "umap", group.by = "sample")

Markers <- FindAllMarkers(MM, min.pct = 0.25, only.pos = TRUE)
write.csv(Markers, "MM_Markers.csv", row.names = T)

# 
FeaturePlot(MM, features = c("Ccr2", "Vcan", "Cd14", "Lyz2", "Fn1", "Itgam"))
FeaturePlot(MM, features = c("Csf1r", "Adgre1", "Mertk", "Apoe", "Lpl"))
FeaturePlot(MM, features = c("Mpo", "Elane", "Prtn3", "Ctsg"))

# Markers <- FindMarkers(obj, assay = "RNA", ident.1 = "Intermediate mono-macro")

# Highlight the MM 1 cells in obj
highlight <- WhichCells(MM, expression = cellType == "Unannotated")
highlight <- WhichCells(MM, idents = 2)
DimPlot(obj, reduction = "umap", cells.highlight = highlight)
# 应该把6和7去掉，因为他们的方向很奇怪

# Cd68是区分粒细胞和monocyte的主要marker，Csf1r在mono/macro都高表达。
FeaturePlot(obj, reduction = "umap", features = c("Csf1r", "Ly6c2", "Cd68"))
VlnPlot(obj, features = c("Csf1r", "Ly6c2", "Cd68"), pt.size = 0)

# 说明2群不是Cycling promyelocyte
FeaturePlot(obj, reduction = "umap", features = c("Mpo", "Elane", "Prtn3"))
FeaturePlot(obj, features = c("Ccr2", "Vcan", "Cd14", "Lyz2", "Fn1", "Itgam"))

MM$cellType <- "Unannotated"
MM$cellType[MM$seurat_clusters %in% c("0")]  <- "Inflammatory Monocyte"
MM$cellType[MM$seurat_clusters %in% c("1")]  <- "Promyelocyte"
MM$cellType[MM$seurat_clusters %in% c("2")]  <- "Cycling Monocyte"
MM$cellType[MM$seurat_clusters %in% c("3")]  <- "cDC2" # MHC-II Hihg APC
MM$cellType[MM$seurat_clusters %in% c("4")]  <- "Monocyte-derived Macrophages"
MM$cellType[MM$seurat_clusters %in% c("5")]  <- "Cycling_MKI67hi"
# MM$cellType[MM$seurat_clusters %in% c("6")]  <- "Neutrophil_Mature_CampNgp"
# MM$cellType[MM$seurat_clusters %in% c("7")]  <- "Neutrophil_Cxcr2_Retnlg"
MM$cellType[MM$seurat_clusters %in% c("8")]  <- "cDC1"
# MM$cellType[MM$seurat_clusters %in% c("9")]  <- "B_cell"
MM$cellType[MM$seurat_clusters %in% c("10")] <- "IFN-stimulated myeloid"
# MM$cellType[MM$seurat_clusters %in% c("11")] <- "Macrophage_Retnla_Fcrls"

DimPlot(MM, reduction = "umap", group.by = "cellType")

removed <- WhichCells(MM, expression = cellType == "Unannotated")

obj@meta.data[colnames(MM), "SubType"] <- MM@meta.data[colnames(MM), "cellType"]

DimPlot(obj, reduction = "umap", group.by = "SubType")

saveRDS(MM, "Base/Mono.rds")

Gran <- readRDS("Base/Gran.rds")
obj@meta.data[colnames(Gran), "SubType"] <- Gran@meta.data[colnames(Gran), "cellType"]
rm(Gran); gc()

obj$SubType[is.na(obj$SubType)] <- obj$MainType[is.na(obj$SubType)]
DimPlot(obj, reduction = "umap", group.by = "SubType")

DimPlot(obj, reduction = "umap", 
        cells.highlight = WhichCells(obj, expression = SubType == "Cycling_MKI67hi"))

# 将Cycling_MKI67hi和Unannotated从obj中移除
obj <- subset(obj, subset = !SubType %in% c("Cycling_MKI67hi", "Unannotated", "Granulocytes"))

DimPlot(obj, reduction = "umap", group.by = "SubType")
obj$group <- gsub("\\d+$", "", obj$sample)

obj$SubType[obj$SubType == "resident-like Macrophage"] <- "Resident-like Macrophage"

DimPlot(obj, reduction = "umap", group.by = "MainType")
DimPlot(obj, reduction = "umap", group.by = "SubType")
DimPlot(obj, reduction = "umap", 
        cells.highlight = WhichCells(obj, expression = SubType == "Granulocytes"))


obj$MainType[obj$SubType == "Resident-like Macrophage"] <- "Macrophage"
obj$MainType[obj$MainType == "Intermediate mono-macro"] <- "Monocyte"
obj$MainType[obj$MainType == "Erythrocytes"] <- "Erythrocyte"
obj$MainType[obj$MainType == "Granulocytes"] <- "Granulocyte"

obj$MainType[obj$MainType == "Intermediate mono-macro"] <- "Monocyte"
obj$MainType[obj$MainType == "Erythrocytes"] <- "Erythrocyte"
obj$MainType[obj$MainType == "Granulocytes"] <- "Granulocyte"
obj$MainType[obj$SubType == "Promyelocyte"] <- "Granulocyte"


pal_main <- c(
  "B cell"              = "#1F77B4",  # 蓝
  "Plasma Cell"         = "#D62728",  # 红
  "Mast Cell Precursor" = "#E377C2",  # 粉
  "Erythrocyte"         = "#8C1D40",  # 酒红（避免和Plasma撞）
  "Granulocyte"         = "#FFD166",  # 浅金色：有颜色但视觉不刺眼
  "Monocyte"            = "#2CA02C",  # 绿
  "Macrophage"          = "#0B7D62",  # 深蓝绿（与Monocyte同系但更深）
  "Dendritic Cell"      = "#FF7F0E",  # 橙
  "NK cell"             = "#9467BD",  # 紫
  "NK-like T cell"      = "#17BECF"  # 青（和 NK 区分）
  )

obj$MainType <- factor(obj$MainType, levels = names(pal_main))

DimPlot(
  obj, reduction = "umap",
  group.by = "MainType",
  cols = pal_main,
  raster = TRUE,
  label.size = 2
) + ggtitle("")


# 堆叠柱状图
proportion_df <- obj@meta.data %>% 
  select(sample, MainType) %>% 
  dplyr::group_by(sample, MainType) %>% 
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(proportion = count / sum(count))

ggplot(proportion_df, aes(x = sample, y = proportion, fill = MainType)) +
  geom_col(position = "stack", width = 0.7, col = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = pal_main, drop = FALSE) +
  theme_classic() +
  labs(x = "sample", y = "Proportion", fill = "MainType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
rm(proportion_df)


proportion_df <- obj@meta.data %>% 
  select(group, MainType) %>% 
  dplyr::group_by(group, MainType) %>% 
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  dplyr::group_by(group) %>% 
  dplyr::mutate(proportion = count / sum(count))

ggplot(proportion_df, aes(x = group, y = proportion, fill = MainType)) +
  geom_col(position = "stack", width = 0.7, col = "black") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = pal_main, drop = FALSE) +
  theme_classic() +
  labs(x = "group", y = "Proportion", fill = "MainType") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
rm(proportion_df)

# DotPlot
Bcell_marker <- c("Cd79a", "Cd79b", "Cd19", "Ms4a1", "Cd22", "Bank1", "Fcer2a")
plasma_marker <- c("Jchain")
mast_marker <- c("Cpa3", "Gata2")
ery_marker <- c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Gata1", "Klf1")
gran_marker <- c("Ly6g", "S100a8", "S100a9", "Lcn2", "Camp", "Ltf", "Cd177")
mono_marker <- c("Ly6c2","Ccr2")
macro_marker <- c("Csf1r", "Cd68", "Adgre1")
DC_marker <- c("Flt3", "Tcf4", "Siglech", "Bst2")
NK_marker <- c("Ncr1", "Klrk1", "Klrb1c", "Nkg7", "Prf1")
Tcell_marker <- c("Cd3e", "Cd8a", "Cd8b1")

DotPlot(obj, features = c(Bcell_marker, plasma_marker, mast_marker,
                          ery_marker, gran_marker, mono_marker, 
                          macro_marker, DC_marker, NK_marker,
                          Tcell_marker), group.by = "MainType") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    panel.grid.major = element_line(color = "grey90")
  ) + scale_color_gradientn(colors = viridis::viridis(10))

DimPlot(obj, reduction = "umap", group.by = "SubType")

# B cell, Gran, Mono的子群体区分。

