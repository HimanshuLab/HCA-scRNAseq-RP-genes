celltype <- "CD8T"

#Input of tissue seurat objects and subsetting the cell type
seurat <- readRDS("Lymph_06_04.rds")
Tissue_name <- "Lymph"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Lymph <- subset(seurat, subset = celltype == "CD8T")
Celltype_Lymph$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "CD8T")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Spleen_06_04.rds")
Tissue_name <- "spleen"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_spleen <- subset(seurat, subset = celltype == "CD8T")
Celltype_spleen$orig.ident <- paste0(Tissue_name,"_",celltype)

merged_seurat <- merge(
  x = Celltype_Marrow,
  y = list(Celltype_spleen,Celltype_Lymph),
  add.cell.ids = c("Marrow","Spleen","Lymph")
)

table(merged_seurat@meta.data$orig.ident)

file_name <- paste0(celltype,"_merged.rds")
saveRDS(merged_seurat,file_name)

#checking the correctness of celltype using markers
Idents(merged_seurat) <- merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(merged_seurat,features = c("CD3D","CD3E","CD8A","CD8B","NKG7"))
plot1



