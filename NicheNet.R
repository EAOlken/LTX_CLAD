#nichenet
```{r}

library(nichenetr)
library(Seurat)
library(tidyverse)

cell_counts <- table(LTx_Seurat2$unk, LTx_Seurat2$group)
keep_celltypes <- rownames(cell_counts)[apply(cell_counts, 1, function(x) all(x >= 3))]

LTx_Seurat2 <- subset(LTx_Seurat2, subset = unk %in% keep_celltypes)


ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))


Idents(LTx_Seurat2) <- LTx_Seurat2$unk
DefaultAssay(LTx_Seurat2) <- "RNA"

condition_oi <-  "BOS"
#condition_oi <-  "RAS"
condition_reference <- "Stable"

seurat_obj_receiver <- subset(LTx_Seurat2, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

vis_ligand_aupr <- ligand_activities %>%
    filter(test_ligand %in% best_upstream_ligands) %>%
    dplyr::select(test_ligand, aupr_corrected) %>%          # explicit select
    arrange(aupr_corrected) %>%
    column_to_rownames("test_ligand") %>%
    as.matrix()

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkblue") + 
        theme(axis.text.x.top = element_blank()))  

DotPlot(subset(LTx_Seurat2, unk %in% sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") + 
    coord_flip() +
    scale_y_discrete(position = "right")+ theme(axis.text.x = element_text(angle = 90))


(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
   theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))


library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(nichenetr)   # if not installed, install from GitHub: devtools::install_github("saeyslab/nichenetr")
DefaultAssay(LTx_Seurat2) <- "RNA"
Idents(LTx_Seurat2) <- LTx_Seurat2$unk
table(LTx_Seurat2$group)

# Do receiver & BOS/Stable exist?
receiver <- "CD38+ ATII"  
"CD38+ ATII" %in% unique(Idents(LTx_Seurat2))
c("BOS","Stable") %in% unique(LTx_Seurat2$group)

#receiver <- "Activated"
#receiver <- "Pro-fibrotic"
# Count per receiver x group
table(Idents(LTx_Seurat2), LTx_Seurat2$group)[receiver, ]

# Keep only sender cell types with >=3 cells per group
keep_celltypes <- rownames(cell_counts)[apply(cell_counts[, c("BOS", "Stable")], 1, function(x) all(x >= 3))]
# Subset your Seurat object to those
LTx_filtered <- subset(LTx_Seurat2, subset = unk %in% keep_celltypes)
Idents(LTx_filtered) <- LTx_filtered$unk

table(Idents(LTx_filtered), LTx_filtered$group)

# Run wrapper (will fail if receiver has <3 cells in a condition)
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = LTx_filtered,
  receiver = receiver,
  condition_colname = "group",
  condition_oi = "RAS",
  condition_reference = "Stable",
  assay_oi = "RNA",
  weighted_networks = weighted_networks,
  lr_network = lr_network,
  ligand_target_matrix = ligand_target_matrix
)
```
