setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS/RDS_files")
celltype <- "CD8T_NKG7"

# seurat <- readRDS("Bladder_16_04.rds")
# Tissue_name <- "Bladder"
# levels(seurat)
# seurat@meta.data$celltype <- Idents(seurat)
# Celltype_bladder <- subset(seurat, subset = celltype == "Macrophage")
# Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Bile_15_04.rds")
levels(seurat)
Tissue_name <- "Bile"
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bile <- subset(seurat, subset = celltype == "CD8T_NKG7+")
Celltype_bile$orig.ident <- paste0(Tissue_name,"_",celltype)

# seurat <- readRDS("esophagous_15_04.rds")
# Tissue_name <- "oesophagus"
# levels(seurat)
# seurat@meta.data$celltype <- Idents(seurat)
# Celltype_eso <- subset(seurat, subset = celltype == "Macrophage")
# Celltype_eso$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Liver_12_04.rds")
Tissue_name <- "Liver"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_liver <- subset(seurat, subset = celltype == "CD8+NKT")
Celltype_liver$orig.ident <- paste0(Tissue_name,"_",celltype)

# seurat <- readRDS("Muscle_05_04_recent.rds")
# Tissue_name <- "Muscle"
# levels(seurat)
# seurat@meta.data$celltype <- Idents(seurat)
# Celltype_Muscle <- subset(seurat, subset = celltype == "Endothelial")
# Celltype_Muscle$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "CD8T")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)

# seurat <- readRDS("Skin_15_04.rds")
# Tissue_name <- "Skin"
# levels(seurat)
# seurat@meta.data$celltype <- Idents(seurat)
# Celltype_Skin <- subset(seurat, subset = celltype == "Endothelial")
# Celltype_Skin$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "CD8T")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "CD8T_NKG7+")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Blood_19_04.rds")
Tissue_name <- "Blood"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Blood <- subset(seurat, subset = celltype == "CD8+ NKT")
Celltype_Blood$orig.ident <- paste0(Tissue_name,"_",celltype)

# seurat <- readRDS("Lymph_06_04.rds")
# Tissue_name <- "Lymph"
# levels(seurat)
# seurat@meta.data$celltype <- Idents(seurat)
# Celltype_Lymph <- subset(seurat, subset = celltype == "Endothelial")
# Celltype_Lymph$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "CD8T NKG7+")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("SmallIntestine_11_06.rds")
Tissue_name <- "SI"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_SI <- subset(seurat, subset = celltype == "CD8T")
Celltype_SI$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Spleen_06_04.rds")
Tissue_name <- "spleen"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_spleen <- subset(seurat, subset = celltype == "CD8+ NKT")
Celltype_spleen$orig.ident <- paste0(Tissue_name,"_",celltype)

#merged_seurat <- merge(x = Celltype_Trachea, y = c(Celltype_Stomach, Celltype_SI, Celltype_Skin, Celltype_Rectum, Celltype_Muscle, Celltype_liver, Celltype_bile, Celltype_bladder))



merged_seurat <- merge(
  x = Celltype_Marrow,
  y = list(Celltype_Rectum, Celltype_liver, Celltype_Blood, Celltype_bile,Celltype_SI,Celltype_spleen,Celltype_Stomach,Celltype_Trachea),
  add.cell.ids = c("Marrow", "Rectum", "Liver", "Blood", "Bile","SI","Spleen", "Stomach","Trachea")
)

table(merged_seurat@meta.data$orig.ident)

VlnPlot(merged_seurat,features = c("CD8A","NKG7"))

#, Celltype_Muscle,Celltype_Skin, Celltype_Stomach, Celltype_spleen, Celltype_SI, Celltype_Marrow, Celltype_Lymph, Celltype_Blood,
file_name <- paste0(celltype,"_merged.rds")
saveRDS(merged_seurat,file_name)

merged_seurat <- readRDS("Endothelial_merged.rds")
Idents(merged_seurat) <- merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(merged_seurat,features = c("CD3D","CD3E","CD8A","CD8B","NKG7"))
plot1



