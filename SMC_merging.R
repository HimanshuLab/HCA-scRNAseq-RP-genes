setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS")
celltype <- "SMC"

seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "Smoothmuscle")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Bile_15_04.rds")
levels(seurat)
Tissue_name <- "Bile"
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bile <- subset(seurat, subset = celltype == "Smooth muscle")
Celltype_bile$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("esophagous_15_04.rds")
Tissue_name <- "oesophagus"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_eso <- subset(seurat, subset = celltype == "Smooth muscle")
Celltype_eso$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Liver_12_04.rds")
Tissue_name <- "Liver"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_liver <- subset(seurat, subset = celltype == "Smoothmuscle")
Celltype_liver$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Muscle_05_04_recent.rds")
Tissue_name <- "Muscle"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Muscle <- subset(seurat, subset = celltype == "Smooth muscle")
Celltype_Muscle$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "SMC")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Skin_15_04.rds")
Tissue_name <- "Skin"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Skin <- subset(seurat, subset = celltype == "Smooth muscle")
Celltype_Skin$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "Smoothmuscle")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)


merged_seurat <- merge(
  x = Celltype_Rectum,
  y = list(Celltype_liver, Celltype_Muscle, Celltype_bile,Celltype_bladder,Celltype_eso,Celltype_Skin,Celltype_Trachea),
  add.cell.ids = c("Rectum", "Liver", "Muscle", "Bile","Bladder","Oesophagus","Skin","Trachea")
)

table(merged_seurat@meta.data$orig.ident)

new_merged_seurat <- merge(
  x = Celltype_Rectum,
  y = list(Celltype_Muscle, Celltype_bladder,Celltype_eso,Celltype_Skin,Celltype_Trachea),
  add.cell.ids = c("Rectum", "Muscle", "Bladder","Oesophagus","Skin","Trachea")
)


#, Celltype_Muscle,Celltype_Skin, Celltype_Stomach, Celltype_spleen, Celltype_SI, Celltype_Marrow, Celltype_Lymph, Celltype_Blood,
file_name <- paste0(celltype,"_merged.rds")
saveRDS(new_merged_seurat,file_name)


Idents(new_merged_seurat) <- new_merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(new_merged_seurat,features = c("ACTA2","MYL9","TAGLN","MYH11"))
plot1



