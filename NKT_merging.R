setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS")
celltype <- "NKT"

seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "NKT")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Bile_15_04.rds")
levels(seurat)
Tissue_name <- "Bile"
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bile <- subset(seurat, subset = celltype == "NKT")
Celltype_bile$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("esophagous_15_04.rds")
Tissue_name <- "oesophagus"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_eso <- subset(seurat, subset = celltype == "NKT")
Celltype_eso$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Liver_12_04.rds")
Tissue_name <- "Liver"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_liver <- subset(seurat, subset = celltype == "NKT")
Celltype_liver$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Muscle_05_04_recent.rds")
Tissue_name <- "Muscle"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Muscle <- subset(seurat, subset = celltype == "NKT")
Celltype_Muscle$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Skin_15_04.rds")
Tissue_name <- "Skin"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Skin <- subset(seurat, subset = celltype == "NKT")
Celltype_Skin$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "NKT")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "NKT")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "NKT")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Spleen_06_04.rds")
Tissue_name <- "spleen"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_spleen <- subset(seurat, subset = celltype == "NKT")
Celltype_spleen$orig.ident <- paste0(Tissue_name,"_",celltype)

merged_seurat <- merge(
  x = Celltype_Marrow,
  y = list(Celltype_liver, Celltype_bile,Celltype_spleen,Celltype_Stomach,Celltype_Trachea, Celltype_bladder, Celltype_eso,Celltype_Muscle, Celltype_Skin),
  add.cell.ids = c("Marrow", "Liver", "Bile","Spleen", "Stomach","Trachea","Bladder","Oesophagus","Muscle","Skin")
)

table(merged_seurat@meta.data$orig.ident)

new_merged_seurat <- merge(
  x = Celltype_Marrow,
  y = list(Celltype_liver, Celltype_bile,Celltype_spleen,Celltype_Trachea, Celltype_bladder, Celltype_eso,Celltype_Muscle),
  add.cell.ids = c("Marrow", "Liver", "Bile","Spleen","Trachea","Bladder","Oesophagus","Muscle")
)

#, Celltype_Muscle,Celltype_Skin, Celltype_Stomach, Celltype_spleen, Celltype_SI, Celltype_Marrow, Celltype_Lymph, Celltype_Blood,
file_name <- paste0(celltype,"_merged.rds")
saveRDS(new_merged_seurat,file_name)

Idents(new_merged_seurat) <- new_merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(new_merged_seurat,features = c("CD3D","CD3E","CD8A","CD8B","NKG7"))
plot1



