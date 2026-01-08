setwd("~/Projects/ImmuOmics_XL")

library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(CellChat)
library(circlize)

# configuration
samples <- unique(obj$sample)
for (sample_use in samples) {
  print(sample_use)
  
  celltype_col <- "MainType"   
  species      <- "mouse" 
  outdir       <- file.path("V01/CellChat_single", sample_use)
  dir.create(outdir, recursive = TRUE, showWarnings = TRUE)
  
  # 单样本 + 用SCT@data
  obj1 <- subset(obj, subset = sample == sample_use)
  data.input <- GetAssayData(obj1, assay = "SCT", layer = "data")
  meta <- data.frame(labels = obj1@meta.data[[celltype_col]], row.names = colnames(data.input))
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat@DB <- if (species == "mouse") CellChatDB.mouse else CellChatDB.human
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat) # most costy
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, file.path(outdir, paste0("cellchat_", sample_use, ".rds")))
  gc()
}

rm(data.input, mat, mat2, meta, obj1, celltype_col, i, outdir, sample_use, samples, species); gc()

# visualization
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")

pathway_show <- if ("ANNEXIN" %in% cellchat@netP$pathways) "ANNEXIN" else cellchat@netP$pathways[1]
netVisual_aggregate(cellchat, signaling = pathway_show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathway_show)
netVisual_bubble(cellchat, signaling = pathway_show, remove.isolate = FALSE)

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") +
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")



########### 测试可用 ########
library(igraph)
library(ggraph)
library(ggplot2)

mat <- cellchat@net$count
mat[is.na(mat)] <- 0
diag(mat) <- 0

gs <- table(cellchat@idents)
gs <- gs[rownames(mat)]  # 对齐顺序

edges <- as.data.frame(as.table(mat))
edges <- edges[edges$Freq > 0, ]
colnames(edges) <- c("from","to","weight")

nodes <- data.frame(name = rownames(mat), size = as.numeric(gs))

g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

ggraph(g, layout = "circle") +
  geom_edge_fan(aes(width = weight, alpha = weight), show.legend = FALSE) +
  geom_node_point(aes(size = size), show.legend = FALSE) +
  geom_node_text(aes(label = name), vjust = 1.2, size = 3) +
  scale_edge_width(range = c(0.2, 3)) +
  theme_void() +
  ggtitle("Number of interactions")

rm(mat, gs, edges, nodes, g); gc()



###########
pathways.show <- c("CXCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


groupSize <- table(cellchat@idents)  # 建议保留名字
netVisual_circle2(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                       label.edge = FALSE, title.name = "Weight")

rm(groupSize, pathway.show, ptm); gc()


##### Merge analysis test #####
AD1_cc <- readRDS("CellChat_single/AD1/cellchat_AD1.rds")
AD3_cc <- readRDS("CellChat_single/AD3/cellchat_AD3.rds")
AS1_cc <- readRDS("CellChat_single/AS1/cellchat_AS1.rds")
Control2_cc <- readRDS("CellChat_single/Control2/cellchat_Control2.rds")

cellchat_merged <- mergeCellChat(list(AS1 = AS1_cc, Control2 = Control2_cc),
                                 add.names = c("AS1","Control2"))

# cellchat_merged <- mergeCellChat(list(AD1 = AD1_cc, Control2 = Control2_cc),
#                                  add.names = c("AD1","Control2"))

# cellchat_merged <- mergeCellChat(list(Control2 = Control2_cc, AD1 = AD1_cc),
#                                  add.names = c("Control2","AD1"))

# p1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))
compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "weight")
# rm(p1, p2)


# netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)

netVisual_heatmap(cellchat_merged, comparison = c(1,2), measure = "weight")



#### AD1 analysis ####
sample <- "AD3"
cellchat <- readRDS(paste0("CellChat_single/", sample, "/cellchat_", sample, ".rds"))


groupSize <- table(cellchat@idents)  # 建议保留名字
netVisual_circle2(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                  label.edge = FALSE, title.name = "Weight")
netVisual_circle2(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                  label.edge = FALSE, title.name = "Count")

rm(groupSize, pathway.show, ptm); gc()

cellchat@netP$pathways

# 通路
pathways.show <- c("CCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.ICAM <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.ICAM[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 7, targets.use = seq_len(nlevels(cellchat@idents)), 
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = seq_len(nlevels(cellchat@idents)), targets.use = 7, 
                 remove.isolate = FALSE)

ptm = Sys.time()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred inte
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_scatter(cellchat)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht2

VlnPlot(obj, features = c("Thy1"), group.by = "MainType", assay = "SCT")


# 通路模式分析
library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "incoming") # 非常费时间
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

netAnalysis_dot(cellchat, pattern = "incoming")
