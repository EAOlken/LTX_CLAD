
```{r}
library(CellChat)
library(Seurat)
library(patchwork)


obj.list <- SplitObject(LTx_Seurat2, split.by = "group")
stablechat <- obj.list$Stable
boschat <- obj.list$BOS
rasschat <- obj.list$RAS

stablechat$newhlca <- droplevels(stablechat$newhlca)
levels(stablechat$newhlca)
boschat$newhlca <- droplevels(boschat$newhlca)
levels(boschat$newhlca)
rasschat$newhlca <- droplevels(rasschat$newhlca)
levels(rasschat$newhlca)

data_inputstable <- GetAssayData(stablechat, assay = "RNA", slot = "data")  # normalized data
data_inputbos <- GetAssayData(boschat, assay = "RNA", slot = "data")  # normalized data
data_inputras <- GetAssayData(rasschat, assay = "RNA", slot = "data")  # normalized data

cellchatstable <- createCellChat(object = data_inputstable, meta = stablechat@meta.data, group.by = "newhlca")
cellchatbos <- createCellChat(object = data_inputbos, meta = boschat@meta.data, group.by = "newhlca")
cellchatras <- createCellChat(object = data_inputras, meta = rasschat@meta.data, group.by = "newhlca")

CellChatDB <- CellChatDB.human
cellchatstable@DB <- CellChatDB
cellchatbos@DB    <- CellChatDB
cellchatras@DB    <- CellChatDB

# Preprocessing STABLE
cellchatstable <- subsetData(cellchatstable)  # subset expression data of genes in the database
cellchatstable <- identifyOverExpressedGenes(cellchatstable)
cellchatstable <- identifyOverExpressedInteractions(cellchatstable)
cellchatstable <- computeCommunProb(cellchatstable)
cellchatstable <- filterCommunication(cellchatstable, min.cells = 10)  # adjust threshold
cellchatstable <- computeCommunProbPathway(cellchatstable)
cellchatstable <- aggregateNet(cellchatstable)
#visualization
netVisual_circle(cellchatstable@net$count, vertex.weight = as.numeric(table(stablechat$new)), weight.scale = T, label.edge = F)

cellchatstable <- netAnalysis_computeCentrality(cellchatstable, slot.name = "netP")
#paperplot
netAnalysis_signalingRole_scatter(cellchatstable)

# Preprocessing BOS
cellchatbos <- subsetData(cellchatbos)  # subset expression data of genes in the database
cellchatbos <- identifyOverExpressedGenes(cellchatbos)
cellchatbos <- identifyOverExpressedInteractions(cellchatbos)
cellchatbos <- computeCommunProb(cellchatbos)
cellchatbos <- filterCommunication(cellchatbos, min.cells = 10)  # adjust threshold
cellchatbos <- computeCommunProbPathway(cellchatbos)
cellchatbos <- aggregateNet(cellchatbos)

cellchatbos <- netAnalysis_computeCentrality(cellchatbos, slot.name = "netP")
#netAnalysis_signalingRole_scatter(cellchatbos)



# Preprocessing RAS
cellchatras <- subsetData(cellchatras)  # subset expression data of genes in the database
cellchatras <- identifyOverExpressedGenes(cellchatras)
cellchatras <- identifyOverExpressedInteractions(cellchatras)
cellchatras <- computeCommunProb(cellchatras)
cellchatras <- filterCommunication(cellchatras, min.cells = 10)  # adjust threshold
cellchatras <- computeCommunProbPathway(cellchatras)
cellchatras <- aggregateNet(cellchatras)
cellchatras <- netAnalysis_computeCentrality(cellchatras, slot.name = "netP")
#netAnalysis_signalingRole_scatter(cellchatras)

library(NMF)
library(ggalluvial)

selectK(cellchatbos, pattern = "outgoing")
nPatterns = 4
cellchatbos <- identifyCommunicationPatterns(cellchatbos, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchatbos, pattern = "outgoing")


#merged
object.list <- list(BOS = cellchatbos,RAS = cellchatras)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#interaction strength
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD38+ ATII")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Effector-memory CD8+")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))



cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "BOS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "BOS",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "RAS",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
```
