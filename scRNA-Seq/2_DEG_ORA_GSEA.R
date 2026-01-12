# ==================================
# DGE (differential gene expression)
# ORA (over representing analysis)
# GSEA (gene set enrichment analysis)
# ==================================

# ========== Library ==========
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
obj2DEG <- function(obj, type_cir, type = NULL, group_cir, group01 = NULL, group02 = NULL,
                    write_2_csv = FALSE, output_file = NULL) {
  
  # DESCRIPTION: get the DEG table from a Seurat obj
  
  md <- obj[[]]
  cells_use <- rownames(md)[md[[type_cir]] == type & md[[group_cir]] %in% c(group01, group02)]
  sub_obj <- subset(obj, cells = cells_use)
  
  deg <- FindMarkers(
    sub_obj,
    assay = "SCT",
    group.by = group_cir,
    ident.1 = group01,
    ident.2 = group02,
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    recorrect_umi = FALSE
  )
  
  deg_gsea <- FindMarkers(
    sub_obj,
    assay = "SCT",
    group.by = group_cir,
    ident.1 = group01,
    ident.2 = group02,
    test.use = "wilcox",
    min.pct = 0,
    logfc.threshold = 0,
    return.thresh = 1,
    recorrect_umi = FALSE
  )
  
  deg <- rbind(deg[deg$avg_log2FC > 0, ][order(-deg[deg$avg_log2FC > 0, "avg_log2FC"]), , drop=FALSE],
               deg[deg$avg_log2FC < 0, ][order( deg[deg$avg_log2FC < 0, "avg_log2FC"]), , drop=FALSE])
  
  if (write_2_csv){
    stopifnot(!is.null(output_file))
    write.csv(deg, output_file)
  }
  
  return(list(deg = deg, deg_gsea = deg_gsea))
}

GO_gsea <- function(deg_gsea, OrgDb, ont = "BP", write_2_csv = FALSE, output_file = NULL) {
  
  # DESCRIPTION: do the GSEA in GO database from no filtered DEG table
  
  geneList <- deg_gsea$avg_log2FC
  names(geneList) <- rownames(deg_gsea)
  geneList <- sort(geneList, decreasing = TRUE)
  
  gene_df <- bitr(names(geneList),
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = OrgDb)
  
  geneList2 <- sort(sapply(split(geneList[gene_df$SYMBOL], gene_df$ENTREZID), max), decreasing = TRUE)
  
  gsea_go <- gseGO(
    geneList     = geneList2,
    OrgDb        = OrgDb,
    keyType      = "ENTREZID",
    ont          = ont,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )
  
  gsea_go_df <- as.data.frame(gsea_go)
  
  gsea_go_df <- rbind(
    gsea_go_df[gsea_go_df$NES > 0, ][order(-gsea_go_df[gsea_go_df$NES > 0, "NES"]), , drop=FALSE],
    gsea_go_df[gsea_go_df$NES < 0, ][order( gsea_go_df[gsea_go_df$NES < 0, "NES"]), , drop=FALSE]
  )
  
  if (write_2_csv){
    stopifnot(!is.null(output_file))
    write.csv(gsea_go_df, output_file)
  }
  
  return(list(gsea_go = gsea_go, gsea_go_df = gsea_go_df))
}

GO_ora <- function(deg, OrgDb, ont = "BP", write_2_csv = FALSE, output_file = NULL) {
  
  # DESCRIPTION: do the ORA in GO database from filtered DEG table
  
  gene <- rownames(deg)  
  gene_df <- clusterProfiler::bitr(gene,
                                   fromType = "SYMBOL",
                                   toType   = "ENTREZID",
                                   OrgDb    = OrgDb)
  
  ora_go <- clusterProfiler::enrichGO(
    gene          = unique(gene_df$ENTREZID),
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  ora_go_df <- as.data.frame(ora_go)
  
  if (write_2_csv){
    stopifnot(!is.null(output_file))
    write.csv(ora_go_df, output_file, row.names = FALSE)
  }
  
  return(list(ora_go = ora_go, ora_go_df = ora_go_df))
}

GO_ora_2 <- function(deg, OrgDb, ont = "BP", write_2_csv = FALSE, output_prefix = NULL) {
  
  # DESCRIPTION: do the ORA in GO database from filtered DEG table according to UP and DOWN
  
  padj_cut = 0.05
  lfc_cut = 0.25
  
  up_gene   <- rownames(deg)[deg$p_val_adj < padj_cut & deg$avg_log2FC >  lfc_cut]
  down_gene <- rownames(deg)[deg$p_val_adj < padj_cut & deg$avg_log2FC < -lfc_cut]
  
  up_df <- clusterProfiler::bitr(up_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  dn_df <- clusterProfiler::bitr(down_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)

  ora_up <- clusterProfiler::enrichGO(
    gene          = unique(up_df$ENTREZID),
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  ora_dn <- clusterProfiler::enrichGO(
    gene          = unique(dn_df$ENTREZID),
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  ora_up_df <- as.data.frame(ora_up)
  ora_dn_df <- as.data.frame(ora_dn)
  ora_up_df$Direction <- "Up"
  ora_dn_df$Direction <- "Down"
  
  ora_all_df <- rbind(ora_up_df, ora_dn_df)
  
  ora_all_df$Direction <- factor(ora_all_df$Direction, levels = c("Up", "Down"))
  ora_all_df <- ora_all_df[order(ora_all_df$Direction, -ora_all_df$FoldEnrichment), ]

  if (write_2_csv) {
    stopifnot(!is.null(output_prefix))
    write.csv(ora_up_df,  paste0(output_prefix, "_ORA_UP.csv"),   row.names = FALSE)
    write.csv(ora_dn_df,  paste0(output_prefix, "_ORA_DOWN.csv"), row.names = FALSE)
    write.csv(ora_all_df, paste0(output_prefix, "_ORA_UPDOWN.csv"), row.names = FALSE)
  }
  
  return(list(ora_up = ora_up, ora_dn = ora_dn, ora_all_df = ora_all_df))
}

DEG_volcanoplot <- function(deg, plot_title = NULL, savePNG = FALSE, output_file = NULL) {
  
  # DESCRIPTION: draw the volcano plot of deg table
  
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
  
  top_lab <- df %>%
    filter(sig != "Not sig") %>%
    arrange(p_val_adj) %>%
    head(10)
  
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
    ggtitle(plot_title)
  
  if (savePNG) {
    stopifnot(!is.null(output_file))
    ggsave(filename = output_file,
           plot = p, width = 6, height = 5, dpi = 300)
  }
  
  return(p)
}

enrichment_plot <- function(enrichment, GO_term, savePNG = FALSE, output_file = NULL) {
  p <- gseaplot2(enrichment, 
            geneSetID = GO_term, 
            title = GO_term,
            pvalue_table = TRUE)
  
  if (savePNG) {
    stopifnot(!is.null(output_file))
    ggsave(filename = output_file,
           plot = p, width = 10, height = 6)
  }
  
  return(p)
}

Run_DEG_ORA_GSEA <- function(obj, type_cir, type, group_cir, group01, group02, output_dir) {
  
  # Total pipeline of DEG, ORA and GSEA from a Seurat obj
  
  # type_cir <- "MainType", type <- "Macrophage"
  # group_cir <- "sample", group01 <- "AD1", group02 <- "Control2"
  # output_dir <- "V02"
  
  cat(paste0(">>>>>>>>>>>", group01, "vs", group02, "_", type, ">>>>>>>>>>>>\n"))
  
  cat("start calcuate DEG...\n")
  deg_result <- obj2DEG(obj = obj, type_cir = type_cir, type = type,
                        group_cir = group_cir, group01 = group01, group02 = group02)
  deg <- deg_result$deg
  deg_gsea <- deg_result$deg_gsea
  
  if (is.null(deg) || nrow(deg) == 0) {
    cat("No DEG found. Skip volcano/ORA/GSEA and return.\n")
    return(invisible(NULL))
  }
  
  output_dir <- paste0(output_dir, "/", group01, "vs", group02, "_", type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("start save DEG result...\n")
  file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_filteredDEG.csv")
  cat(paste0("DEG filename : ", file_name, "\n"))
  write.csv(deg, file_name)
  
  cat("start do GSEA...\n")
  gsea_result <- GO_gsea(deg_gsea = deg_gsea, OrgDb = org.Mm.eg.db, ont = "BP")
  gsea <- gsea_result$gsea_go
  gsea_df <- gsea_result$gsea_go_df
  rm(gsea_result)
  
  if (!is.null(gsea_df) && nrow(gsea_df) > 0) {
    cat("start save GSEA result...\n")
    file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_GO_gsea.csv")
    cat(paste0("GSEA filename : ", file_name, "\n"))
    write.csv(gsea_df, file_name)
  } else {
    cat("No GSEA terms -> skip saving GSEA table.\n")
  }
  
  cat("start do ORA...\n")
  # ora_result <- GO_ora(deg = deg, OrgDb = org.Mm.eg.db, ont = "BP")
  ora_result <- GO_ora_2(deg = deg, OrgDb = org.Mm.eg.db, ont = "BP")
  # ora <- ora_result$ora_go
  ora_up <- ora_result$ora_up
  ora_down <- ora_result$ora_dn
  # ora_df <- ora_result$ora_go_df
  ora_df <- ora_result$ora_all_df
  rm(ora_result)
  
  if (!is.null(ora_df) && nrow(ora_df) > 0) {
    cat("start save ORA result...\n")
    file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_GO_ora.csv")
    cat(paste0("ORA filename : ", file_name, "\n"))
    write.csv(ora_df, file_name, row.names = FALSE)
  } else {
    cat("No ORA terms -> skip saving ORA.\n")
  }
  
  if (!is.null(deg) && nrow(deg) > 0 && all(c("avg_log2FC","p_val_adj") %in% colnames(deg))) {
    cat("start draw volcano plot...\n")
    deg_volcano <- DEG_volcanoplot(deg, plot_title = paste0(type, " ", group01, "vs", group02))
    cat("strat save volcano plot...\n")
    file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_DEG_volcano.png")
    cat(paste0("volcano plot file : ", file_name, "\n"))
    ggsave(file_name, plot = deg_volcano, width = 6, height = 5, dpi = 300)
  }
  
  # cat("start draw ora bubble plot...\n")
  # ora_bubble <- dotplot(ora, showCategory = 15)
  # cat("strat save volcano plot...\n")
  # file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_GO_ora_bubble.png")
  # cat(paste0("ora bubble plot file : ", file_name, "\n"))
  # ggsave(file_name, plot = ora_bubble, width = 6, height = 8, dpi = 300)
  
  cat("start draw ora bubble plot (up/down)...\n")
  p1 <- if (is.null(ora_up)   || nrow(as.data.frame(ora_up))   == 0) ggplot() + theme_void() + ggtitle(paste0("Up in ", group01, " (none)"))   else dotplot(ora_up,   showCategory = 15) + ggtitle(paste0("Up in ", group01))
  p2 <- if (is.null(ora_down) || nrow(as.data.frame(ora_down)) == 0) ggplot() + theme_void() + ggtitle(paste0("Down in ", group01, " (none)")) else dotplot(ora_down, showCategory = 15) + ggtitle(paste0("Down in ", group01))
  ora_bubble <- (p1 | p2) + patchwork::plot_layout(guides = "collect")
  cat("start save ora bubble plot...\n")
  file_name <- file.path(output_dir, paste0(group01, "vs", group02, "_", type, "_GO_ora_bubble.png"))
  cat(paste0("ora bubble plot file : ", file_name, "\n"))
  ggplot2::ggsave(file_name, plot = ora_bubble, width = 12, height = 8, dpi = 300)
  
  cat("start draw gsea bubble plot...\n")
  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    gsea_bubble <- dotplot(gsea, x = "NES", color = "qvalue", split = ".sign", showCategory = 15) +
      facet_grid(. ~ .sign)
    gsea_bubble$layers[[1]]$mapping$size <- rlang::sym("setSize")
    gsea_bubble <- gsea_bubble + scale_color_distiller(palette = "RdPu", direction = -1) +
      theme(axis.text.x = element_text(size = 6))
    
    cat("start save gsea bubble plot...\n")
    file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_GO_gsea_bubble.png")
    cat(paste0("gsea bubble plot file : ", file_name, "\n"))
    ggsave(file_name, plot = gsea_bubble, width = 8, height = 12, dpi = 300)
  } else {
    cat("No GSEA terms -> skip gsea bubble plot.\n")
  }
  
  cat("start draw gsea plot...\n")
  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    show_n <- min(5, nrow(as.data.frame(gsea)))
    gsea_plot <- gseaplot2(gsea, 1:show_n, title = "Test", pvalue_table = TRUE)
    
    cat("start save gsea plot...\n")
    file_name <- paste0(output_dir, "/", group01, "vs", group02, "_", type, "_GO_gsea.png")
    cat(paste0("gsea plot file : ", file_name, "\n"))
    ggsave(file_name, plot = gsea_plot, width = 8, height = 5, dpi = 300)
  } else {
    cat("No GSEA terms -> skip gseaplot2.\n")
  }
  
  cat(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
}


# ========== Analysis ==========
# 1. Test: AD1 vs AS1 Macrophage
Run_DEG_ORA_GSEA(obj, type_cir = "MainType", type = "Macrophage", 
                 group_cir = "sample", group01 = "AD1", group02 = "Control2", 
                 output_dir = "./V02")

# 2. Loop to get all the result
for (type in unique(obj$MainType)){
  for (group01 in unique(obj$sample)){
    for (group02 in unique(obj$sample)){
      print(paste0("=============", type, "_", group01, "_", group02, "============"))
      if (group01 != group02){
        Run_DEG_ORA_GSEA(obj, type_cir = "MainType", type = type, 
                         group_cir = "sample", group01 = group01, group02 = group02, 
                         output_dir = "./V02/MainType_sample")
      }
      print("===========================")
    }
  }
}
