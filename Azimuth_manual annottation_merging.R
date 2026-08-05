library(Azimuth)
LTx_Seurat <- RenameIdents(object = LTx_Seurat, "0"="Fibroblast", "1" = "ATI","2"="ATII","3" = "T Cells","4"  = "Endothelial","5"="Macrophages","6"="IM","11"= "Plasma","13" = "SMC","14" = "Basal Cells","15"="B Cells","16"="LEC","17"="Mast Cells")

LTx_Seurat$CellType <- Idents(LTx_Seurat)

LTx_Seurat<-RunAzimuth(LTx_Seurat, "lungref")
Idents(LTx_Seurat) <- LTx_Seurat$predicted.ann_level_1

#correlation map

manual_labels <- LTx_Seurat2$CellType 
azimuth_labels <- LTx_Seurat2$predicted.ann_finest_level

library(ComplexHeatmap)
library(circlize)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Create confusion matrix
confusion_mat <- table(manual_labels,azimuth_labels)

# Convert to percentage by row
confusion_mat_pct <- prop.table(confusion_mat, margin = 1) * 100

# Create labels matrix (same dimensions)
labels <- matrix(
    sprintf("%.1f%%", confusion_mat_pct),
    nrow = nrow(confusion_mat_pct),
    ncol = ncol(confusion_mat_pct),
    dimnames = dimnames(confusion_mat_pct)
)

# Draw heatmap
Heatmap(
    as.matrix(confusion_mat_pct),
    name = "Percentage",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = colorRamp2(c(0, 50, 100), c("white", "skyblue", "darkblue")),
    column_title = "Manual Annotation",
    row_title = "Azimuth Annotation", row_names_gp = gpar(fontsize = 8),        # Row labels
  column_names_gp = gpar(fontsize = 8),
)
