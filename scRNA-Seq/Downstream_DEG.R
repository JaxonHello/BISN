# ======================
# Downstream Analysis
# ======================

# ====== Library ======
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(rlang)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

# ====== Configuration ======
work_dir <- "/Users/jaxonhe/Projects/ImmuOmics_XL/V01"
sc_obj <- "/Users/jaxonhe/Projects/ImmuOmics_XL/V01/sc_obj.RData"

setwd(work_dir)

# ====== Read Data ======
load(sc_obj)
rm(work_dir, sc_obj)

# ====== Downstream Analysis ======

# ====== 2) DGE ======
# ====== 2) MainType
# DEG
obj <- PrepSCTFindMarkers(obj)

maintype <- "NK cell"
sample01 <- "AD1"
sample02 <- "Control2"

for (maintype in unique(obj$MainType)){
  print(maintype)
  for (sample01 in unique(obj$sample)) {
    for (sample02 in unique(obj$sample)) {
      if (sample01 != sample02) {
        print(paste0(sample01, "vs", sample02))
        mac <- subset(obj, subset = MainType == maintype & sample %in% c(sample01, sample02))
        deg <- FindMarkers(
          mac,
          assay = "SCT",
          group.by = "sample",
          ident.1 = sample01,
          ident.2 = sample02,
          test.use = "wilcox",
          min.pct = 0.1,
          logfc.threshold = 0.25,
          recorrect_umi = FALSE
        )
        
        write.csv(deg, paste0("DEG/", sample01, "_", sample02, "_", maintype, ".csv"))
        rm(mac); gc()
        
        # GO enrichment
        geneList <- deg$avg_log2FC
        names(geneList) <- rownames(deg)
        geneList <- sort(geneList, decreasing = TRUE)
        
        gene_df <- bitr(names(geneList),
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org.Mm.eg.db)
        
        geneList2 <- geneList[gene_df$SYMBOL]
        names(geneList2) <- gene_df$ENTREZID
        geneList2 <- sort(geneList2, decreasing = TRUE)
        rm(geneList, gene_df); gc()
        
        gsea_go <- gseGO(
          geneList     = geneList2,
          OrgDb        = org.Mm.eg.db,
          keyType      = "ENTREZID",
          ont          = "BP",
          minGSSize    = 10,
          maxGSSize    = 500,
          pvalueCutoff = 0.05,
          verbose      = FALSE
        )
        
        gsea_go_df <- as.data.frame(gsea_go)
        
        if (nrow(gsea_go_df) > 0){
          # print(dotplot(gsea_go, showCategory = 15, color = "NES"))
          write.csv(gsea_go_df, paste0("DEG/", sample01, "_", sample02, "_", maintype, "_GOBP.csv"))
        } else {
          print(paste0(sample01, "_", sample02, "_", maintype, ":无显著富集"))
        }
        gc()
      }
    }
  }
}


# DEG plot
# volcano plot
fns <- list.files("DEG", pattern="\\.csv$", full.names=TRUE)[!grepl("GOBP\\.csv$", list.files("DEG", pattern="\\.csv$", full.names=TRUE))]
fns

for (fn in fns){
  print(fn)
  deg <- read.csv(fn, row.names = 1)
  
  sample1 <- strsplit(sub("\\.csv$", "", basename(fn)), "_")[[1]][1]
  sample2 <- strsplit(sub("\\.csv$", "", basename(fn)), "_")[[1]][2]
  celltype <- paste(strsplit(sub("\\.csv$", "", basename(fn)), "_")[[1]][-(1:2)], collapse = "_")
  
  lfc_cut <- 0.25
  padj_cut <- 0.05
  
  df <- deg %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      neglog10_padj = -log10(p_val_adj + 1e-300),
      sig = case_when(
        p_val_adj < padj_cut & avg_log2FC >  lfc_cut ~ "Up in Sample1",
        p_val_adj < padj_cut & avg_log2FC < -lfc_cut ~ "Up in Sample2",
        TRUE ~ "Not sig"
      )
    )
  
  deg$gene <- rownames(deg)
  deg$status <- "NS"
  deg$status[deg$p_val_adj < 0.05 & deg$avg_log2FC >  0.25] <- "Up"
  deg$status[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.25] <- "Down"
  
  
  # 标注几个最显著的基因（可改 10/20）
  top_lab <- df %>%
    filter(sig != "Not sig") %>%
    arrange(p_val_adj) %>%
    head(10)
  
  ggplot(df, aes(x = avg_log2FC, y = neglog10_padj)) +
    geom_point(aes(color = sig), size = 1) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = top_lab,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(x = "avg_log2FC (AD1 vs Control2)", y = "-log10(P adj)", color = "") +
    theme_classic()
  
  p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj)) +
    geom_point(aes(color = sig), size = 1) +
    scale_color_manual(values = c("Up in Sample1" = "red", "Up in Sample2" = "blue", "Not sig" = "grey")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = top_lab,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(x = "avg_log2FC", y = "-log10(P adj)", color = "") +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(paste0(sample1, " vs ", sample2, " in ", celltype))
  
  ggsave(filename = paste0("DEG/", sample1, "_vs_", sample2, "_", celltype, "_volcano.png"),
         plot = p, width = 6, height = 5, dpi = 300)
}

# GSEA bubble plot
# GO enrichment
fns <- list.files("DEG", pattern="\\.csv$", full.names=TRUE)[!grepl("GOBP\\.csv$", list.files("DEG", pattern="\\.csv$", full.names=TRUE))]
fns

for (fn in fns){
  print(fn)
  
  sample1 <- strsplit(sub("\\.csv$", "", basename(fn)), "_")[[1]][1]
  sample2 <- strsplit(sub("\\.csv$", "", basename(fn)), "_")[[1]][2]
  celltype <- paste(strsplit(sub("\\.csv$", "", basename(fn)), "_")[[1]][-(1:2)], collapse = "_")
  
  deg <- read.csv(fn, row.names = 1)
  
  geneList <- deg$avg_log2FC
  names(geneList) <- rownames(deg)
  geneList <- sort(geneList, decreasing = TRUE)
  
  gene_df <- bitr(names(geneList),
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Mm.eg.db)
  
  geneList2 <- geneList[gene_df$SYMBOL]
  names(geneList2) <- gene_df$ENTREZID
  geneList2 <- sort(geneList2, decreasing = TRUE)
  rm(geneList, gene_df); gc()
  
  gsea_go <- gseGO(
    geneList     = geneList2,
    OrgDb        = org.Mm.eg.db,
    keyType      = "ENTREZID",
    ont          = "BP",
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )
  
  gsea_go_df <- as.data.frame(gsea_go)
  
  if (nrow(gsea_go_df) > 0){
    # print(dotplot(gsea_go, showCategory = 15, color = "NES"))
    write.csv(gsea_go_df, paste0("DEG/", sample01, "_", sample02, "_", maintype, "_GOBP.csv"))
  } else {
    print(paste0(sample01, "_", sample02, "_", maintype, ":无显著富集"))
  }
  gc()
}




lfc_cut <- 0.25
padj_cut <- 0.05

df <- deg %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    neglog10_padj = -log10(p_val_adj + 1e-300),
    sig = case_when(
      p_val_adj < padj_cut & avg_log2FC >  lfc_cut ~ "Up in Sample1",
      p_val_adj < padj_cut & avg_log2FC < -lfc_cut ~ "Up in Sample2",
      TRUE ~ "Not sig"
    )
  )

deg$gene <- rownames(deg)
deg$status <- "NS"
deg$status[deg$p_val_adj < 0.05 & deg$avg_log2FC >  0.25] <- "Up"
deg$status[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.25] <- "Down"


# 标注几个最显著的基因（可改 10/20）
top_lab <- df %>%
  filter(sig != "Not sig") %>%
  arrange(p_val_adj) %>%
  head(10)

ggplot(df, aes(x = avg_log2FC, y = neglog10_padj)) +
  geom_point(aes(color = sig), size = 1) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_lab,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(x = "avg_log2FC (AD1 vs Control2)", y = "-log10(P adj)", color = "") +
  theme_classic()

p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj)) +
  geom_point(aes(color = sig), size = 1) +
  scale_color_manual(values = c("Up in Sample1" = "red", "Up in Sample2" = "blue", "Not sig" = "grey")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_lab,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(x = "avg_log2FC", y = "-log10(P adj)", color = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample1, " vs ", sample2, " in ", celltype))

ggsave(filename = paste0("DEG/", sample1, "_vs_", sample2, "_", celltype, "_volcano.png"),
       plot = p, width = 6, height = 5, dpi = 300)
