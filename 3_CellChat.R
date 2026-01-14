# ==================================
# Cell Chat analysis
# ==================================

# ========== Library ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(CellChat)
  library(circlize)
  library(igraph)
  library(ggraph)
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
netVisual_circle2 <- function(net, color.use = NULL, title.name = NULL, sources.use = NULL,
                              targets.use = NULL, idents.use = NULL, remove.isolate = FALSE,
                              top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL,
                              vertex.size.max = NULL, vertex.label.cex = 1, vertex.label.color = "black",
                              edge.weight.max = NULL, edge.width.max = 3, alpha.edge = 0.5,
                              label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8,
                              edge.curved = 0.5, shape = "circle", margin = 0.2,
                              vertex.size = NULL, arrow.width = 1, arrow.size = 0.2) {
  
  # deps
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  stopifnot(requireNamespace("ggraph", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("scales", quietly = TRUE))
  stopifnot(requireNamespace("grid", quietly = TRUE))
  stopifnot(requireNamespace("reshape2", quietly = TRUE))
  
  if (!is.null(vertex.size)) warning("'vertex.size' is deprecated. Use `vertex.weight`")
  
  # ---- threshold top edges ----
  net[is.na(net)] <- 0
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  
  # ---- subset source/target/idents ----
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
    if (is.null(rownames(net))) stop("The input weighted matrix should have rownames!")
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- cells.level[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- cells.level[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) idents.use <- cells.level[idents.use]
      df.net <- dplyr::filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  
  # ---- remove isolates ----
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx  <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, , drop = FALSE]
      net <- net[, -idx, drop = FALSE]
    }
  }
  
  # ---- build graph ----
  g <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE, diag = TRUE)
  
  # ---- palette (try to mimic CellChat scPalette) ----
  if (is.null(color.use)) {
    # 先尝试拿 CellChat 的 scPalette（如果你装的是 CellChat）
    if ("CellChat" %in% loadedNamespaces() || requireNamespace("CellChat", quietly = TRUE)) {
      scPal <- try(getFromNamespace("scPalette", "CellChat"), silent = TRUE)
      if (!inherits(scPal, "try-error")) {
        color.use <- scPal(length(igraph::V(g)))
      }
    }
    if (is.null(color.use)) {
      color.use <- grDevices::hcl.colors(length(igraph::V(g)), "Set2")
    }
  }
  
  # ---- vertex weight scaling (match original) ----
  vnames <- igraph::V(g)$name
  
  if (length(vertex.weight) == 1) {
    vertex.weight <- rep(vertex.weight, length(vnames))
    names(vertex.weight) <- vnames
  } else {
    if (is.null(names(vertex.weight))) names(vertex.weight) <- vnames
    vertex.weight <- vertex.weight[vnames]
  }
  
  if (is.null(vertex.size.max)) {
    vertex.size.max <- if (length(unique(vertex.weight)) == 1) 5 else 15
  }
  if (is.null(vertex.weight.max)) vertex.weight.max <- max(vertex.weight, na.rm = TRUE)
  
  vertex.size <- vertex.weight / vertex.weight.max * vertex.size.max + 5
  
  # ---- edges data.frame ----
  edges <- igraph::as_data_frame(g, what = "edges")
  colnames(edges)[colnames(edges) == "from"] <- "source"
  colnames(edges)[colnames(edges) == "to"]   <- "target"
  if (!"weight" %in% colnames(edges)) edges$weight <- 0
  
  # edge width like original
  if (is.null(edge.weight.max)) edge.weight.max <- max(edges$weight, na.rm = TRUE)
  if (edge.weight.max == 0) edge.weight.max <- 1
  
  if (isTRUE(weight.scale)) {
    edges$width <- 0.3 + edges$weight / edge.weight.max * edge.width.max
  } else {
    edges$width <- 0.3 + edge.width.max * edges$weight
  }
  
  # edge color = source node color + alpha
  node_col <- setNames(color.use[seq_along(vnames)], vnames)
  edges$edge_col <- grDevices::adjustcolor(node_col[edges$source], alpha.f = alpha.edge)
  
  # optional edge labels
  if (isTRUE(label.edge)) {
    edges$label <- round(edges$weight, 1)
  } else {
    edges$label <- NA_character_
  }
  
  edges$is_loop <- edges$source == edges$target
  
  # ---- nodes df ----
  nodes <- data.frame(
    name  = vnames,
    size  = vertex.size,
    color = node_col[vnames],
    stringsAsFactors = FALSE
  )
  print(nodes)
  
  # ---- layout circle coords (mimic in_circle + label placement) ----
  n <- nrow(nodes)
  ang <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  dir_deg <- ang * 180 / pi
  nodes$dir_deg <- dir_deg
  x <- cos(ang); y <- sin(ang)
  
  # scale coords similar to scale()
  x <- as.numeric(scale(x))
  y <- as.numeric(scale(y))
  
  # label distance similar to original: vertex.weight/max + 2
  label.dist <- vertex.weight / max(vertex.weight, na.rm = TRUE) + 2
  label.scale <- 1  # 经验系数，让标签离圆更远一点
  lx <- x * label.scale
  ly <- y * label.scale
  
  # label angle similar (clockwise) for readability
  angle_deg <- (ang / (2*pi)) * 360
  angle_deg <- ifelse(angle_deg > 180, angle_deg - 180, angle_deg)  # 翻转到 0-180
  hjust <- ifelse(cos(ang) >= 0, 0, 1)
  
  nodes$x <- x; nodes$y <- y
  nodes$lx <- lx; nodes$ly <- ly
  nodes$angle <- angle_deg
  nodes$hjust <- hjust
  nodes$dir_deg <- (ang * 180 / pi) %% 360 
  
  edges$direction <- nodes$dir_deg[match(edges$source, nodes$name)]
  edges$span <- 80
  
  # ---- build igraph with our node/edge attrs for ggraph ----
  g2 <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  igraph::E(g2)$width   <- edges$width
  igraph::E(g2)$edge_col <- edges$edge_col
  igraph::E(g2)$label   <- edges$label
  
  # ---- plot ----
  p <- ggraph::ggraph(g2, layout = "manual", x = nodes$x, y = nodes$y) +
    ggraph::geom_edge_fan(
      ggplot2::aes(filter = !is_loop, width = I(width), colour = I(edge_col)),
      alpha = 0.6,
      show.legend = FALSE,
      arrow = grid::arrow(length = grid::unit(arrow.size, "cm"), type = "closed")
    ) +
    ggraph::geom_edge_loop(
      ggplot2::aes(filter = is_loop, width = I(width), colour = I(edge_col),
                   direction = direction, span = span),
      alpha = 0.6,
      show.legend = FALSE,
      strength = edge.curved
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, size = size.Freq),
      colour = nodes$color,
      show.legend = FALSE
    ) +
    ggplot2::scale_size(range = c(3, 8), guide = "none") +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = lx, y = ly, label = name, angle = 0, hjust = hjust),
      family = "Helvetica",
      size = 5 * vertex.label.cex,
      colour = vertex.label.color
    ) +
    ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(plot.margin = ggplot2::margin(25,25,25,25)) +
    ggplot2::coord_equal(clip = "off")
  
  if (!is.null(title.name)) p <- p + ggplot2::ggtitle(title.name)
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(p)
}

# ========== Analysis ==========
Run_CellChat_Plot <- function(cellchat, outdir, prefix) {
  
  # DESCRIPTION: run the total cellchat plot pipeline
  
  outdir <- paste0(outdir, "/", prefix)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  groupSize <- table(cellchat@idents) 
  p <- netVisual_circle2(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                         label.edge = FALSE, title.name = "Weight")
  file_name <- paste0(outdir, "/cellchat_weight_circle.png")
  ggsave(filename = file_name, plot = p, width = 8, height = 5, dpi = 300)
  
  p <- netVisual_heatmap(cellchat, measure = "weight")
  file_name <- paste0(outdir, "/cellchat_weight_heatmap.png")
  png(file_name, width = 6, height = 5, units = "in", res = 300)
  ComplexHeatmap::draw(p)
  dev.off()
  
  # single pathway
  for (pathway.show in cellchat@netP$pathways) {
    print(pathway.show)
    
    subdir <- paste0(outdir, "/pathway/", pathway.show)
    if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
    
    # 1) chord
    png(file.path(subdir, paste0(pathway.show, "_chord.png")),
        width = 6, height = 5, units = "in", res = 300)
    netVisual_aggregate(cellchat, signaling = pathway.show, layout = "chord")
    dev.off()
    
    # 2) heatmap
    ht <- netVisual_heatmap(cellchat, signaling = pathway.show, color.heatmap = "Reds")
    png(file.path(subdir, paste0(pathway.show, "_heatmap.png")),
        width = 7, height = 6, units = "in", res = 300)
    ComplexHeatmap::draw(ht)
    dev.off()
    
    # 3) contribution
    p <- netAnalysis_contribution(cellchat, signaling = pathway.show)
    ggsave(file.path(subdir, paste0(pathway.show, "_contribution.png")),
           plot = p, width = 7, height = 5, dpi = 300)
    
    # 4) singaling roles
    png(file.path(subdir, paste0(pathway.show, "_role_network.png")),
        width = 8, height = 2.5, units = "in", res = 300)
    netAnalysis_signalingRole_network(cellchat, signaling = pathway.show,
                                      width = 8, height = 2.5, font.size = 10)
    dev.off()
    
    # 5) singaling roles
    gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathway.show)
    ggplot2::ggsave(file.path(subdir, paste0(pathway.show, "_role_scatter.png")),
                    plot = gg2, width = 6, height = 5, dpi = 300)
  }
  
  # total pathway
  subdir <- paste0(outdir, "/pathway")
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
  
  for (i in seq_len(nlevels(cellchat@idents))) {
    celltype <- levels(cellchat@idents)[i]
    print(celltype)
    
    file_name <- file.path(subdir, paste0(celltype, "_source.png"))
    png(file_name,
        width = 4, height = 6, units = "in", res = 300)
    p <- netVisual_bubble(cellchat, sources.use = i, targets.use = seq_len(nlevels(cellchat@idents)), remove.isolate = FALSE)
    print(p)
    dev.off()
    print(file_name)
    
    
    file_name <- file.path(subdir, paste0(celltype, "_target.png"))
    png(file_name,
        width = 4, height = 6, units = "in", res = 300)
    p <- netVisual_bubble(cellchat, sources.use = seq_len(nlevels(cellchat@idents)), targets.use = i, 
                          remove.isolate = FALSE)
    print(p)
    dev.off()
    print(file_name)
    
  }
  
  # singaling roles total
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  ggplot2::ggsave(file.path(outdir, "signalingRole_scatter.png"), plot = gg1, width = 6, height = 5, dpi = 300)
  
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
  png(file.path(outdir, "signalingRole_out_in.png"), width = 16, height = 6, units = "in", res = 300)
  ComplexHeatmap::draw(ht1 + ht2)
  dev.off()
  
}

for (sample in c("AD1", "AD3", "AS1", "AS2", "Aged1", "Control2")) {
  print(sample)
  cellchat <- readRDS(paste0("CellChat_single/", sample, "/cellchat_", sample, ".rds"))
  Run_CellChat_Plot(cellchat = cellchat, outdir = "V02/3_CellChat", prefix = sample)
}


Run_CellChat_Merged_Plot <- function(cellchat01, cellchat02, names, outdir, prefix) {
  
  # DESCRIPTION: run the cellchat comparison pipeline
  
  outdir <- paste0(outdir, "/", prefix)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  cellchat_merged <- mergeCellChat(list(cellchat01, cellchat02),
                                   add.names = names)
  
  ptm = Sys.time()
  gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
  p <- gg1 + gg2
  ggplot2::ggsave(file.path(outdir, "interaction_comparison.png"), plot = p, width = 10, height = 5, dpi = 300)
  
  
  gg1 <- netVisual_heatmap(cellchat_merged)
  gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
  ht <- gg1 + gg2
  png(file.path(outdir, "interaction_comparison_heatmap.png"),
      width = 10, height = 6, units = "in", res = 300)
  ComplexHeatmap::draw(ht)
  dev.off()
  
  for (ct in levels(cellchat_merged@meta$labels)) {
    
    print(ct)
    subdir <- paste0(outdir, "/", ct)
    if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
    
    p <- netAnalysis_signalingChanges_scatter(cellchat_merged,
                                         idents.use = ct,
                                         signaling.exclude = "MIF")
    ggplot2::ggsave(file.path(subdir, "singaling_change.png"), plot = p, width = 6, height = 5, dpi = 300)
    
    
    
    file_name <- file.path(subdir, paste0(ct, "_source.png"))
    png(file_name,
        width = 10, height = 6, units = "in", res = 300)
    gg1 <- netVisual_bubble(cellchat_merged,
                            sources.use = ct, targets.use = levels(cellchat_merged@idents),
                            comparison = c(1,2), max.dataset = 2,
                            title.name = "Increased signaling in dataset2",
                            angle.x = 45, remove.isolate = TRUE
    )
    
    gg2 <- netVisual_bubble(cellchat_merged,
                            sources.use = ct, targets.use = levels(cellchat_merged@idents),
                            comparison = c(1,2), max.dataset = 1,
                            title.name = "Decreased signaling in dataset2",
                            angle.x = 45, remove.isolate = TRUE
    )
    print(gg1 + gg2)
    dev.off()
    
    file_name <- file.path(subdir, paste0(ct, "_target.png"))
    png(file_name,
        width = 10, height = 6, units = "in", res = 300)
    gg1 <- netVisual_bubble(cellchat_merged,
                            sources.use = levels(cellchat_merged@idents), targets.use = ct,
                            comparison = c(1,2), max.dataset = 2,
                            title.name = "Increased signaling in dataset2",
                            angle.x = 45, remove.isolate = TRUE
    )
    
    gg2 <- netVisual_bubble(cellchat_merged,
                            sources.use = levels(cellchat_merged@idents), targets.use = ct,
                            comparison = c(1,2), max.dataset = 1,
                            title.name = "Decreased signaling in dataset2",
                            angle.x = 45, remove.isolate = TRUE
    )
    print(gg1 + gg2)
    dev.off()
  }
  
}

Run_CellChat_Merged_Plot_02 <- function(cellchat01, cellchat02, names, outdir, prefix) {
  
  # ---- helpers ----
  safe_run <- function(desc, expr) {
    tryCatch(
      withCallingHandlers(
        expr,
        warning = function(w) {
          cat("WARN  |", desc, ":", conditionMessage(w), "\n")
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        cat("ERROR |", desc, ":", conditionMessage(e), "\n")
        NULL
      }
    )
  }
  
  safe_png <- function(file, width = 10, height = 6, units = "in", res = 300) {
    ok <- TRUE
    safe_run(paste0("png() open -> ", file), {
      grDevices::png(file, width = width, height = height, units = units, res = res)
    })
    # 如果打开失败，当前设备仍是 null device(1)
    if (grDevices::dev.cur() == 1) ok <- FALSE
    ok
  }
  
  safe_dev_off <- function() {
    # 只有当前不是 null device 才关
    if (grDevices::dev.cur() != 1) {
      safe_run("dev.off()", { grDevices::dev.off() })
    }
  }
  
  # ---- outdir ----
  outdir <- file.path(outdir, prefix)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- merge ----
  cellchat_merged <- safe_run("mergeCellChat()", {
    CellChat::mergeCellChat(list(cellchat01, cellchat02), add.names = names)
  })
  if (is.null(cellchat_merged)) return(invisible(NULL))
  
  # 取细胞类型 levels（你之前 levels(cellchat_merged@idents) 是 NULL，所以这里用 meta$labels）
  ct_levels <- safe_run("get celltype levels", {
    labs <- cellchat_merged@meta$labels
    if (is.null(labs)) character(0) else levels(factor(labs))
  })
  if (is.null(ct_levels) || length(ct_levels) == 0) {
    cat("WARN  | no cell types found in cellchat_merged@meta$labels\n")
    return(invisible(NULL))
  }
  
  # ---- compareInteractions plot ----
  p_cmp <- safe_run("compareInteractions()", {
    gg1 <- CellChat::compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))
    gg2 <- CellChat::compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "weight")
    gg1 + gg2
  })
  safe_run("save interaction_comparison.png", {
    if (!is.null(p_cmp)) {
      ggplot2::ggsave(file.path(outdir, "interaction_comparison.png"),
                      plot = p_cmp, width = 10, height = 5, dpi = 300)
    }
  })
  
  # ---- heatmap compare ----
  safe_run("netVisual_heatmap() + save", {
    ht1 <- CellChat::netVisual_heatmap(cellchat_merged)
    ht2 <- CellChat::netVisual_heatmap(cellchat_merged, measure = "weight")
    ht  <- ht1 + ht2
    
    fn <- file.path(outdir, "interaction_comparison_heatmap.png")
    if (safe_png(fn, width = 10, height = 6)) {
      ComplexHeatmap::draw(ht)
      safe_dev_off()
    }
  })
  
  # ---- per cell type ----
  for (ct in ct_levels) {
    cat("== celltype:", ct, "==\n")
    
    subdir <- file.path(outdir, ct)
    if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    
    # signaling change scatter
    p_sc <- safe_run(paste0("netAnalysis_signalingChanges_scatter(", ct, ")"), {
      CellChat::netAnalysis_signalingChanges_scatter(cellchat_merged,
                                                     idents.use = ct,
                                                     signaling.exclude = "MIF")
    })
    safe_run(paste0("save signaling_change.png (", ct, ")"), {
      if (!is.null(p_sc)) {
        ggplot2::ggsave(file.path(subdir, "singaling_change.png"),
                        plot = p_sc, width = 6, height = 5, dpi = 300)
      }
    })
    
    # bubble: ct as source (increase/decrease)
    safe_run(paste0("netVisual_bubble SOURCE (", ct, ") + save"), {
      fn <- file.path(subdir, paste0(ct, "_source.png"))
      if (safe_png(fn, width = 10, height = 6)) {
        
        gg1 <- safe_run(paste0("bubble source increased (", ct, ")"), {
          CellChat::netVisual_bubble(cellchat_merged,
                                     sources.use = ct,
                                     targets.use = ct_levels,
                                     comparison = c(1,2),
                                     max.dataset = 2,
                                     title.name = "Increased signaling in dataset2",
                                     angle.x = 45,
                                     remove.isolate = TRUE)
        })
        
        gg2 <- safe_run(paste0("bubble source decreased (", ct, ")"), {
          CellChat::netVisual_bubble(cellchat_merged,
                                     sources.use = ct,
                                     targets.use = ct_levels,
                                     comparison = c(1,2),
                                     max.dataset = 1,
                                     title.name = "Decreased signaling in dataset2",
                                     angle.x = 45,
                                     remove.isolate = TRUE)
        })
        
        if (!is.null(gg1) && !is.null(gg2)) print(gg1 + gg2)
        safe_dev_off()
      }
    })
    
    # bubble: ct as target (increase/decrease)
    safe_run(paste0("netVisual_bubble TARGET (", ct, ") + save"), {
      fn <- file.path(subdir, paste0(ct, "_target.png"))
      if (safe_png(fn, width = 10, height = 6)) {
        
        gg1 <- safe_run(paste0("bubble target increased (", ct, ")"), {
          CellChat::netVisual_bubble(cellchat_merged,
                                     sources.use = ct_levels,
                                     targets.use = ct,
                                     comparison = c(1,2),
                                     max.dataset = 2,
                                     title.name = "Increased signaling in dataset2",
                                     angle.x = 45,
                                     remove.isolate = TRUE)
        })
        
        gg2 <- safe_run(paste0("bubble target decreased (", ct, ")"), {
          CellChat::netVisual_bubble(cellchat_merged,
                                     sources.use = ct_levels,
                                     targets.use = ct,
                                     comparison = c(1,2),
                                     max.dataset = 1,
                                     title.name = "Decreased signaling in dataset2",
                                     angle.x = 45,
                                     remove.isolate = TRUE)
        })
        
        if (!is.null(gg1) && !is.null(gg2)) print(gg1 + gg2)
        safe_dev_off()
      }
    })
    
    gc()
  }
  
  invisible(cellchat_merged)
}

Run_CellChat_Merged_Plot_02(cellchat01 = AD1_cc, cellchat02 = AS1_cc, names = c("AD1", "AS1"), 
                            outdir = "V02/3_CellChat", prefix = paste0(c1, "vs", c2))

for (c1 in c("AD1", "AD3")) {
  for (c2 in c("AS1", "AS2", "Aged1", "Control2")) {
    print(paste0(c1, "vs", c2))
    cellchat01 <- readRDS(paste0("CellChat_single/", c1, "/cellchat_", c1, ".rds"))
    cellchat02 <- readRDS(paste0("CellChat_single/", c2, "/cellchat_", c2, ".rds"))
    Run_CellChat_Merged_Plot_02(cellchat01 = cellchat01, cellchat02 = cellchat02, names = c(c1, c2), 
                             outdir = "V02/3_CellChat", prefix = paste0(c1, "vs", c2))
  }
}

for (c1 in c("AS1", "AS2")) {
  for (c2 in c("AD1", "Aged1", "Control2")) {
    print(paste0(c1, "vs", c2))
    cellchat01 <- readRDS(paste0("CellChat_single/", c1, "/cellchat_", c1, ".rds"))
    cellchat02 <- readRDS(paste0("CellChat_single/", c2, "/cellchat_", c2, ".rds"))
    Run_CellChat_Merged_Plot_02(cellchat01 = cellchat01, cellchat02 = cellchat02, names = c(c1, c2), 
                                outdir = "V02/3_CellChat", prefix = paste0(c1, "vs", c2))
  }
}


for (c1 in c("Aged1")) {
  for (c2 in c("Control2")) {
    print(paste0(c1, "vs", c2))
    cellchat01 <- readRDS(paste0("CellChat_single/", c1, "/cellchat_", c1, ".rds"))
    cellchat02 <- readRDS(paste0("CellChat_single/", c2, "/cellchat_", c2, ".rds"))
    Run_CellChat_Merged_Plot_02(cellchat01 = cellchat01, cellchat02 = cellchat02, names = c(c1, c2), 
                                outdir = "V02/3_CellChat", prefix = paste0(c1, "vs", c2))
  }
}

