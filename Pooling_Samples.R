library(Seurat)
library(readxl)
library(dplyr)

# All batches and BOS folders
batches <- paste0("batch", 1:4)
bos_dirs <- paste0("BOS", 1:3)

# -------------------------------
# 2. Loop through all folders
# -------------------------------
seurat_list <- list()

for (batch in batches) {
  for (bos in bos_dirs) {
    
    parent_dir <- file.path(batch, "cellranger", bos, "outs", "per_sample_outs")
    
    # Get all sample IDs inside BOS folder
    if (dir.exists(parent_dir)) {
      sample_ids <- list.dirs(parent_dir, full.names = FALSE, recursive = FALSE)
      
      for (sid in sample_ids) {
        sample_path <- file.path(parent_dir, sid, "count", "sample_filtered_feature_bc_matrix")
        
        if (dir.exists(sample_path)) {
          message("Reading: ", batch, " / ", bos, " / ", sid)
          
          data <- Read10X(sample_path)
          obj <- CreateSeuratObject(counts = data, project = paste(batch, bos, sid, sep = "_"))
          
          # Add metadata
          obj$batch <- batch
          obj$BOS <- bos
          obj$sample <- sid
          
          seurat_list[[paste(batch, bos, sid, sep = "_")]] <- obj
        }
      }
    }
  }
}

# -------------------------------
# 3. Merge needed Seurat objects
# -------------------------------

Metadata <- read_excel("Metadata.xlsx")
bos <- subset(Metadata, Group == c("LTX-BOS"))
ras <- subset(Metadata, Group == c("LTX-RAS"))
stable <- subset(Metadata, Group == c("LTX-stable"))
donor <- subset(Metadata, Group == c("Donor"))

comb <- rbind(stable, donor, bos, ras)

names(seurat_list) <- sapply(strsplit(names(seurat_list), "_"), tail, 1)

seurat_list2 <- seurat_list[names(seurat_list) %in% comb$`Sample ID`]
bos_samples     <- seurat_list2[names(seurat_list2) %in% comb$`Sample ID`[comb$Group == "LTX-BOS"]]
ras_samples     <- seurat_list2[names(seurat_list2) %in% comb$`Sample ID`[comb$Group == "LTX-RAS"]]
stable_samples  <- seurat_list2[names(seurat_list2) %in% comb$`Sample ID`[comb$Group == "LTX-stable"]]
donor_samples  <- seurat_list2[names(seurat_list2) %in% comb$`Sample ID`[comb$Group == "Donor"]]

# Merge inside each group (if more than 1 sample)
if (length(bos_samples) > 1) {
  bos_merged <- merge(bos_samples[[1]], y = bos_samples[-1], add.cell.ids = names(bos_samples), project = "LTX-BOS")
} else {
  bos_merged <- bos_samples[[1]]
}

if (length(ras_samples) > 1) {
 ras_merged <- merge(ras_samples[[1]], y = ras_samples[-1], add.cell.ids = names(ras_samples), project = "LTX-RAS")
} else {
  ras_merged <- ras_samples[[1]]
}

if (length(stable_samples) > 1) {
  stable_merged <- merge(stable_samples[[1]], y = stable_samples[-1], add.cell.ids = names(stable_samples), project = "LTX-stable")
} else {
  stable_merged <- stable_samples[[1]]
}

if (length(donor_samples) > 1) {
  donor_samples <- merge(donor_samples[[1]], y = donor_samples[-1], add.cell.ids = names(donor_samples), project = "Donor")
} else {
  donor_merged <- donor_samples[[1]]
}

stable_merged$group <- "Stable"
donor_merged$group <- "Donor"
ras_merged$group <- "RAS"
bos_merged$group <- "BOS"

stable_donor_bos_ras <- merge(x=stable_merged, y=c(donor_samples,bos_merged,ras_merged))

