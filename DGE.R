# pseudo-bulk 组件DE

# =============================
# 1) Pseudobulk counts by (MainType, sample)
# =============================
DefaultAssay(obj) <- "RNA"
pb_counts <- AggregateExpression(
  obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("MainType", "sample"),
  return.seurat = FALSE
)[["RNA"]]


# ======================================
# 1) Cell composition <-> composition
# ======================================
get_prop_mat <- function(obj, type_col, sample_col, min_cells_per_type = 30) {
  md <- obj@meta.data %>%
    transmute(
      sample = .data[[sample_col]],
      type   = .data[[type_col]]
    ) %>%
    filter(!is.na(sample), !is.na(type))
  
  # 计数：每个 sample 里每个 type 的细胞数
  tab <- md %>%
    count(sample, type, name = "n")
  
  # 过滤太稀有的 type（可选，防止相关性不稳定）
  keep_type <- tab %>%
    group_by(type) %>%
    summarise(total = sum(n), .groups = "drop") %>%
    filter(total >= min_cells_per_type) %>%
    pull(type)
  
  tab <- tab %>% filter(type %in% keep_type)
  
  # 变宽：sample × type
  wide <- tab %>%
    group_by(sample) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(sample, type, prop) %>%
    pivot_wider(names_from = type, values_from = prop, values_fill = 0)
  
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$sample  # 行=sample，列=type
  mat
}
# CLR 变换（适用于组成数据）
clr_transform <- function(mat, pseudocount = 1e-6) {
  mat2 <- mat + pseudocount
  gm <- exp(rowMeans(log(mat2)))       # 每个 sample 的几何均值
  log(mat2 / gm)                       # CLR
}
# 计算相关矩阵 + p 值矩阵（Hmisc::rcorr 很方便）
get_cor_and_p <- function(mat, method = "spearman") {
  suppressPackageStartupMessages(library(Hmisc))
  rc <- Hmisc::rcorr(as.matrix(mat), type = ifelse(method == "pearson","pearson","spearman"))
  list(r = rc$r, p = rc$P)
}
plot_cor_heatmap <- function(r_mat, p_mat = NULL, title = "") {
  df <- as.data.frame(as.table(r_mat))
  colnames(df) <- c("type1", "type2", "r")
  
  if (!is.null(p_mat)) {
    df$p <- as.vector(p_mat)
    df$star <- cut(df$p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***","**","*",""))
  } else {
    df$star <- ""
  }
  
  ggplot(df, aes(type1, type2, fill = r)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = star), size = 3) +
    coord_equal() +
    scale_fill_gradient2(limits = c(-1, 1), midpoint = 0) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(title = title, fill = "corr")
}

mat_main <- get_prop_mat(obj, type_col = "MainType", sample_col = "sample")

# A) 直接用比例做 Spearman
res_main_spear <- get_cor_and_p(mat_main, method = "spearman")
p1 <- plot_cor_heatmap(res_main_spear$r, res_main_spear$p,
                       title = "MainType composition correlation (Spearman on proportions)")

# B) CLR 后做 Pearson（推荐）
mat_main_clr <- clr_transform(mat_main)
res_main_clr <- get_cor_and_p(mat_main_clr, method = "pearson")
p2 <- plot_cor_heatmap(res_main_clr$r, res_main_clr$p,
                       title = "MainType composition correlation (Pearson on CLR)")


df_wide <- data.frame(table(obj$sample, obj$MainType)) %>%
  rename(sample = Var1, MainType = Var2, n = Freq) %>%
  pivot_wider(names_from = MainType, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("sample") %>%
  mutate(across(everything(), ~ .x / sum(.x)))
  
cor_mat <- cor(as.matrix(df_wide), method = "spearman", use = "pairwise.complete.obs")

df_cor <- as.data.frame(as.table(cor_mat)) %>%
  rename(Type1 = Var1, Type2 = Var2, r = Freq)

ggplot(df_cor, aes(Type1, Type2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, low = "#2C7BB6", mid = "white", high = "#D7191C") +
  coord_equal() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9)
  ) +
  labs(fill = "Spearman r")


cor_clr <- data.frame(table(obj$sample, obj$MainType)) %>%
  dplyr::rename(sample = Var1, MainType = Var2, n = Freq) %>%
  tidyr::pivot_wider(names_from = MainType, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("sample") %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / sum(.x))) %>%     # proportion
  as.matrix() %>%
  (\(m) {                                                                    # CLR
    m <- m + 1e-6
    log(m / exp(rowMeans(log(m))))
  })() %>%
  cor(method = "pearson", use = "pairwise.complete.obs")


as.data.frame(as.table(cor_clr)) %>%
  setNames(c("Type1","Type2","r")) %>%
  ggplot(aes(Type1, Type2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0,
                       low = "#2C7BB6", mid = "white", high = "#D7191C") +
  coord_equal() +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.text.y = element_text(size = 9)) +
  labs(fill = "CLR Pearson r")
