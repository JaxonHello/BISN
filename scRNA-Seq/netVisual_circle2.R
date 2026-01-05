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


