celltype <- "Fibroblast"

#Input of tissue seurat objects and subsetting the cell type
seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "Fibroblast")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Bile_15_04.rds")
levels(seurat)
Tissue_name <- "Bile"
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bile <- subset(seurat, subset = celltype == "Fibroblast")
Celltype_bile$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("esophagous_15_04.rds")
Tissue_name <- "oesophagus"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_eso <- subset(seurat, subset = celltype == c("Fibroblast C7+","Fibroblast PTGDS+","Fibroblast PLA2G2A+","Fibroblast APOD+"))
Celltype_eso$orig.ident <- paste0(Tissue_name,"_",celltype)


seurat <- readRDS("Muscle_05_04_recent.rds")
Tissue_name <- "Muscle"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Muscle <- subset(seurat, subset = celltype == "Fibroblast")
Celltype_Muscle$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "Fibroblasts")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Skin_15_04.rds")
Tissue_name <- "Skin"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Skin <- subset(seurat, subset = celltype == "Fibroblast")
Celltype_Skin$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "Fibroblast")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "Fibroblast")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)



merged_seurat <- merge(
  x = Celltype_bladder,
  y = list(Celltype_eso, Celltype_Muscle, Celltype_Rectum,Celltype_Skin,Celltype_Stomach,Celltype_Trachea),
  add.cell.ids = c("Bladder", "Oesophagus", "Muscle", "Rectum","Skin", "Stomach","Trachea")
)

table(merged_seurat@meta.data$orig.ident)

file_name <- paste0(celltype,"_merged.rds")
saveRDS(merged_seurat,file_name)

#checking the correctness of celltype using markers
Idents(merged_seurat) <- merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(merged_seurat,features = c("MMP2","DCN"))
plot1



