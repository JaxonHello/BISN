# ==============================================================================
# Bulk Gene Expression Difference Analysis
# ==============================================================================

# ===== Library ===== 
library(data.table)
library(DESeq2)
library(glmnet)
library(ggplot2)
library(dplyr)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)

# ===== Work dir =====
# Type work dir here
setwd("~/Projects/PhD_Projects/Ding Infor/LQ_1224/TCGA-LUAD.star_counts and clinical")

# ===== Read data =====
# Counts matrix
counts_data <- fread("TCGA-LUAD.star_counts.tsv", data.table = F, check.names = FALSE)
rownames(counts_data) <- counts_data$Ensembl_ID
counts_data <- counts_data[,-1] 

# Phenotype
pheno_data <- fread("TCGA-LUAD.clinical.tsv", data.table = F, check.names = FALSE)
rownames(pheno_data) <- pheno_data$sample

# ===== Preprocess =====
# 1. return to raw counts
counts_data <- 2^counts_data - 1
counts_data <- round(counts_data)
counts_data[counts_data < 0] <- 0

# 2. intersect table
common_samples <- intersect(colnames(counts_data), rownames(pheno_data))
counts_data <- counts_data[, common_samples]
pheno_data <- pheno_data[common_samples, ]
rm(common_samples); gc() # save memory

# 3. filter the "Tumor" and "Control"
target_samples <- which(pheno_data$`sample_type.samples` %in% c("Primary Tumor", "Solid Tissue Normal"))
counts_data <- counts_data[, target_samples]
pheno_data <- pheno_data[target_samples, ]
rm(target_samples); gc()

# 4. identify group
pheno_data$group <- factor(
  ifelse(pheno_data$`sample_type.samples` == "Primary Tumor", "Tumor", "Normal"),
  levels = c("Normal", "Tumor") 
)

# ===== DGE: differential Gene expression =====
# 1. DESeq DGE
colData <- data.frame(row.names = rownames(pheno_data), condition = pheno_data$group)
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
sig_genes <- subset(res, padj < 0.01 & abs(log2FoldChange) > 2)
sig_gene_ids <- rownames(sig_genes)


# ===== Downstream & Visualization =====
# 1. PCA 
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# P1: PCA reduction
p1 <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha=0.8) +
  scale_color_manual(values=c("#377EB8", "#E41A1C")) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  theme_bw() +
  ggtitle("PCA Analysis: LUAD vs Normal")
ggsave("Figure1_PCA.png", p1, width = 6, height = 5)

# 2. DGE
res_df <- as.data.frame(res)
res_df$diff <- "NO"
res_df$diff[res_df$log2FoldChange > 2 & res_df$padj < 0.01] <- "UP"
res_df$diff[res_df$log2FoldChange < -2 & res_df$padj < 0.01] <- "DOWN"

# P2: DGE volcano plot
p2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diff)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_classic() +
  labs(title="Volcano Plot", x="log2(Fold Change)", y="-log10(FDR)")
ggsave("Figure2_Volcano.png", p2, width = 6, height = 5)

# 3. GO enrichment
gene_list <- rownames(sig_genes)
gene_list_clean <- sub("\\..*", "", gene_list)

# P3: Bubble plot
ego <- enrichGO(gene          = gene_list_clean,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)

p3 <- dotplot(ego, showCategory=15, title="GO Enrichment Analysis (Biological Process)")
ggsave("Figure3_GO_Enrichment.png", p3, width = 8, height = 7)

# 4. LASSO Logistic regression
group <- ifelse(pheno_data$`sample_type.samples` == "Primary Tumor", "Tumor", "Normal")
group <- factor(group, levels = c("Normal", "Tumor"))
lasso_data <- assay(vsd)[sig_gene_ids, , drop = FALSE]

x <- t(lasso_data)
x <- as.matrix(x)
storage.mode(x) <- "double"

# 关键：确保 glmnet 能拿到变量名（否则就用 1,2,3…/173,243… 这种编号）
colnames(x) <- rownames(lasso_data)
y <- as.numeric(group) - 1 

set.seed(123) 
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10,
                    standardize = TRUE, type.measure = "deviance")

png("Figure4_LASSO_CV.png", width = 1200, height = 1000, res=150)
plot(cv_fit)
dev.off()

fit <- glmnet(x, y, family = "binomial", alpha = 1, standardize = TRUE)
png("Figure5_LASSO_Path.png", width = 1200, height = 1000, res=150)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE)
dev.off()

lambda_min <- cv_fit$lambda.min
coef_min <- coef(cv_fit, s = "lambda.min")
active_index <- which(coef_min != 0)
active_genes_raw <- rownames(coef_min)[active_index]
active_genes <- active_genes_raw[active_genes_raw != "(Intercept)"]

final_genes_clean <- sub("\\..*", "", active_genes)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = final_genes_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

result_table <- data.frame(Ensembl_ID = active_genes, Gene_Symbol = gene_symbols)

write.csv(result_table, "Result_Lasso_Genes_with_Names.csv")

rm(gsea_go_df, plot_df, celltype, file, gobp_files, group1, group2, p, pathways.show)



