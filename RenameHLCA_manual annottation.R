LTx_Seurat$unk <- as.character(Idents(LTx_Seurat))
LTx_Seurat$unk[Cells(Epi)] <- paste(Idents(Epi))
LTx_Seurat$unk[Cells(Bas)] <- paste(Idents(Bas))
LTx_Seurat$unk[Cells(Fib)] <- paste(Idents(Fib))
LTx_Seurat$unk[Cells(Mac)] <- paste(Idents(Mac))
LTx_Seurat$unk[Cells(TC)] <- paste(Idents(TC))
LTx_Seurat$unk[Cells(cd8)] <- paste(Idents(cd8))
LTx_Seurat$unk[Cells(End)] <- paste(Idents(End))

Idents(LTx_Seurat2) <- LTx_Seurat2$unk
LTx_Seurat2$hlca <- LTx_Seurat2$predicted.ann_finest_level
ltx <- subset(LTx_Seurat2, idents=c("HHIP+ ATII","ATI","Stressed","MUC5B+ Goblet","LEC","Pericyte","EC_Venous","CD4+","EC_aCAP","CD69+ CD8+","Stem-like Basal","Club","NK Cell","TRM_1","Mast Cells","MUC5AC+ Goblet","EC_Immune","COL1A1+ ATII","TRM_2","Proliferating","Alveolar","Lipofibroblast","Stem-like CD8+","TRM CD8+","CD38+ ATII","TREG","Pro-fibrotic","Effector-memory CD8+","Activated","Basal Cells","Aberrant Basoloid","TREM2+ IM","Adventitial","AT2 Progenitor","Primed Basal","Suprabasal","CXCL13+CD4+","TWEAK+TGFB+","AGTR2+ ATII"))

LTx_Seurat2$hlca[Cells(ltx)] <- paste(Idents(ltx))
Idents(LTx_Seurat2) <- LTx_Seurat2$hlca
LTx_Seurat2<- RenameIdents(LTx_Seurat2, "NK cells"="NK Cell", "Alveolar"="Alveolar Fibroblasts","Adventitial"="Adventitial fibroblasts","Activated"="Activated Fibroblasts","Stressed"="Stressed Fibroblasts","Pro-fibrotic"="Pro-fibrotic Fibroblasts", "EC general capillary"="EC_gCAP", "EC aerocyte capillary"="EC_aCAP", "EC venous pulmonary"="EC_Venous","AT1"="ATI")

DimPlot_scCustom(LTx_Seurat2, repel = T, label.size=3, split.by = "group",split_seurat = TRUE, shuffle = FALSE,pt.alpha=0.2)+NoLegend()
