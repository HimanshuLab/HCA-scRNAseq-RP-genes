#Input of tissue seurat objects and subsetting the cell type
celltype <- "T"

seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "T")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Bile_15_04.rds")
levels(seurat)
Tissue_name <- "Bile"
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bile <- subset(seurat, subset = celltype == "T")
Celltype_bile$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Liver_12_04.rds")
Tissue_name <- "Liver"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_liver <- subset(seurat, subset = celltype == "T")
Celltype_liver$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Muscle_05_04_recent.rds")
Tissue_name <- "Muscle"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Muscle <- subset(seurat, subset = celltype == "T cell")
Celltype_Muscle$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "T")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)


seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "T")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "T cell")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Blood_19_04_undefrem_recent.rds")
Tissue_name <- "Blood"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Blood <- subset(seurat, subset = celltype == "T")
Celltype_Blood$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Lymph_06_04.rds")
Tissue_name <- "Lymph"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Lymph <- subset(seurat, subset = celltype == "T")
Celltype_Lymph$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "T")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Small_intetsine_06_04.rds")
Tissue_name <- "SI"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_SI <- subset(seurat, subset = celltype == "T")
Celltype_SI$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Spleen_06_04.rds")
Tissue_name <- "spleen"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_spleen <- subset(seurat, subset = celltype == "T")
Celltype_spleen$orig.ident <- paste0(Tissue_name,"_",celltype)


merged_seurat <- merge(
  x = Celltype_bile,
  y = list(Celltype_Rectum, Celltype_liver, Celltype_Blood, Celltype_Marrow,Celltype_SI,Celltype_spleen,Celltype_Stomach,Celltype_Trachea,Celltype_Lymph, Celltype_Muscle ,Celltype_bladder),
  add.cell.ids = c("Bile","Rectum", "Liver", "Blood", "Marrow","SI","Spleen", "Stomach","Trachea","Lymph","Muscle","Bladder")
)

table(merged_seurat@meta.data$orig.ident)

#tissues with less than 100 cells not considered
new_merged_seurat <- merge(
  x = Celltype_bile,
  y = list(Celltype_Rectum, Celltype_liver, Celltype_Blood, Celltype_Marrow,Celltype_spleen,Celltype_Stomach,Celltype_Trachea,Celltype_Lymph ,Celltype_bladder),
  add.cell.ids = c("Bile","Rectum", "Liver", "Blood", "Marrow","Spleen", "Stomach","Trachea","Lymph","Bladder")
)

file_name <- paste0(celltype,"_merged.rds")
saveRDS(new_merged_seurat,file_name)

#checking the correctness of cell type
Idents(new_merged_seurat) <- new_merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(new_merged_seurat,features = c("CD3D","CD3E","CD8A","CD8B","NKG7"))
plot1



