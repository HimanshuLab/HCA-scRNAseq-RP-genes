celltype <- "Epithelial"

#Input of tissue seurat objects and subsetting the cell type
seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "Epithelial")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("esophagous_15_04.rds")
Tissue_name <- "oesophagus"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_eso <- subset(seurat, subset = celltype == c("Squamous epithelium KRT4+","Squamous epithelium","Squamous epithelium KRT14+","Highpro epithelium","Columnar epithelium"))
Celltype_eso$orig.ident <- paste0(Tissue_name,"_",celltype)


seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "Epithelial")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == c("Basal_Epithelial","Secretory_Epithelial"))
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)


merged_seurat <- merge(
  x = Celltype_eso,
  y = list(Celltype_Rectum,Celltype_bladder,Celltype_Trachea),
  add.cell.ids = c("Oesophagus", "Rectum","Bladder","Trachea")
)

table(merged_seurat@meta.data$orig.ident)
file_name <- paste0(celltype,"_merged.rds")
saveRDS(merged_seurat,file_name)

#checking the correctness of celltype using markers
Idents(merged_seurat) <- merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(merged_seurat,features = c("KRT18","KRT19"))
plot1



