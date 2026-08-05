library(Seurat)
library(SeuratExtend)
library(readxl)

#Scoreing
```{r}
Sens_Apo <- read_excel("Sens_Apo.xlsx") #attached

LTx_Seurat2 <- AddModuleScore(
    object = LTx_Seurat2,
    features = list(Sens_Apo$Apoptosis_GO),
    name = "Apoptosis_Score"
)

VlnPlot2( LTx_Seurat2,
    features = "Apoptosis_Score1",
    group.by = "CellType",
    split.by = "group",
    pt.size = 0,
    assay = "RNA", pt = FALSE,
    stat.method = "wilcox.test"
)
```
