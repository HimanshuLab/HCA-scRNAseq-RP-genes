celltype <- "Plasma"

#Input of tissue seurat objects and subsetting the cell type
seurat <- readRDS("esophagous_15_04.rds")
Tissue_name <- "oesophagus"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_eso <- subset(seurat, subset = celltype == "Plasma")
Celltype_eso$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "Plasma")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "Plasma")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "Plasma")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Blood_19_04_undefrem_recent.rds")
Tissue_name <- "Blood"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Blood <- subset(seurat, subset = celltype == "Plasma")
Celltype_Blood$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "Plasma")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Small_intetsine_06_04.rds")
Tissue_name <- "SI"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_SI <- subset(seurat, subset = celltype == "Plasma")
Celltype_SI$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Spleen_06_04.rds")
Tissue_name <- "spleen"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_spleen <- subset(seurat, subset = celltype == "Plasma")
Celltype_spleen$orig.ident <- paste0(Tissue_name,"_",celltype)

#tissues with less than 100 cells not considered
merged_seurat <- merge(
  x = Celltype_Stomach,
  y = list(Celltype_Rectum, Celltype_spleen),
  add.cell.ids = c("Stomach", "Rectum", "Spleen")
)

table(merged_seurat@meta.data$orig.ident)
file_name <- paste0(celltype,"_merged.rds")
saveRDS(merged_seurat,file_name)

#checking the correctness of celltype using markers
merged_seurat <- readRDS("Endothelial_merged.rds")
Idents(merged_seurat) <- merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(merged_seurat,features = c("CD79A","JCHAIN"))
plot1



