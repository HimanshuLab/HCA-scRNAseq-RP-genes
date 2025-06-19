setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS")
celltype <- "Monocyte"

seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "Monocyte")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Bile_15_04.rds")
levels(seurat)
Tissue_name <- "Bile"
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bile <- subset(seurat, subset = celltype == "Monocyte")
Celltype_bile$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Liver_12_04.rds")
Tissue_name <- "Liver"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_liver <- subset(seurat, subset = celltype == "Monocyte")
Celltype_liver$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "Monocyte")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Blood_19_04_undefrem_recent.rds")
Tissue_name <- "Blood"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Blood <- subset(seurat, subset = celltype == "Monocyte")
Celltype_Blood$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Lymph_06_04.rds")
Tissue_name <- "Lymph"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Lymph <- subset(seurat, subset = celltype == "Monocyte")
Celltype_Lymph$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "Monocytes")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

merged_seurat <- merge(
  x = Celltype_Marrow,
  y = list(Celltype_Stomach, Celltype_liver, Celltype_Blood, Celltype_bile, Celltype_bladder, Celltype_Lymph),
  add.cell.ids = c("Marrow", "Stomach", "Liver", "Blood", "Bile","Bladder","Lymph")
)

table(merged_seurat@meta.data$orig.ident)

new_merged_seurat <- merge(
  x = Celltype_liver,
  y = list(Celltype_bile, Celltype_bladder),
  add.cell.ids = c("Liver","Bile","Bladder")
)


#, Celltype_Muscle,Celltype_Skin, Celltype_Stomach, Celltype_spleen, Celltype_SI, Celltype_Marrow, Celltype_Lymph, Celltype_Blood,
file_name <- paste0(celltype,"_merged.rds")
saveRDS(new_merged_seurat,file_name)

Idents(new_merged_seurat) <- new_merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(new_merged_seurat,features = c("C1QA","S100A8","S100A9","S100A12","VCAN"))
plot1



