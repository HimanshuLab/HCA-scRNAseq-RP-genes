setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS")
celltype <- "Tgd"

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "Tgd")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "Tgd TRDV1+")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Small_intetsine_06_04.rds")
Tissue_name <- "SI"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_SI <- subset(seurat, subset = celltype == "Tgd")
Celltype_SI$orig.ident <- paste0(Tissue_name,"_",celltype)



merged_seurat <- merge(
  x = Celltype_Marrow,
  y = list(Celltype_Rectum, Celltype_SI),
  add.cell.ids = c("Marrow", "Rectum", "SI")
)

table(merged_seurat@meta.data$orig.ident)



#, Celltype_Muscle,Celltype_Skin, Celltype_Stomach, Celltype_spleen, Celltype_SI, Celltype_Marrow, Celltype_Lymph, Celltype_Blood,
file_name <- paste0(celltype,"_merged.rds")
saveRDS(merged_seurat,file_name)

Idents(merged_seurat) <- merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(merged_seurat,features = c("CD3D","CD3E","TRDV1","TRDV2"))
plot1



