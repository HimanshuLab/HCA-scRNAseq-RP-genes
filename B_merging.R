setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS")
celltype <- "B"

seurat <- readRDS("Bladder_16_04.rds")
Tissue_name <- "Bladder"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_bladder <- subset(seurat, subset = celltype == "B")
Celltype_bladder$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Liver_12_04.rds")
Tissue_name <- "Liver"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_liver <- subset(seurat, subset = celltype == "B")
Celltype_liver$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Rectum_06_04.rds")
Tissue_name <- "Rectum"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Rectum <- subset(seurat, subset = celltype == "B")
Celltype_Rectum$orig.ident <- paste0(Tissue_name,"_",celltype)


seurat <- readRDS("Stomach_05_04.rds")
Tissue_name <- "Stomach"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Stomach <- subset(seurat, subset = celltype == "B")
Celltype_Stomach$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Trachea_14_04.rds")
Tissue_name <- "Trachea"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Trachea <- subset(seurat, subset = celltype == "B")
Celltype_Trachea$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Blood_19_04_undefrem_recent.rds")

Tissue_name <- "Blood"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Blood <- subset(seurat, subset = celltype == "B")
Celltype_Blood$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Lymph_06_04.rds")
Tissue_name <- "Lymph"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Lymph <- subset(seurat, subset = celltype == "B")
Celltype_Lymph$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Marrow_17_04.rds")
Tissue_name <- "Marrow"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_Marrow <- subset(seurat, subset = celltype == "B")
Celltype_Marrow$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Small_intetsine_06_04.rds")
Tissue_name <- "SI"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_SI <- subset(seurat, subset = celltype == "B")
Celltype_SI$orig.ident <- paste0(Tissue_name,"_",celltype)

seurat <- readRDS("Spleen_06_04.rds")
Tissue_name <- "spleen"
levels(seurat)
seurat@meta.data$celltype <- Idents(seurat)
Celltype_spleen <- subset(seurat, subset = celltype == "B")
Celltype_spleen$orig.ident <- paste0(Tissue_name,"_",celltype)

merged_seurat <- merge(
  x = Celltype_Rectum,
  y = list(Celltype_liver, Celltype_Blood, Celltype_Marrow,Celltype_SI,Celltype_spleen,Celltype_Stomach,Celltype_Trachea,Celltype_Lymph,Celltype_bladder),
  add.cell.ids = c("Rectum", "Liver", "Blood", "Marrow","SI","Spleen", "Stomach","Trachea","Lymph","Bladder")
)

table(merged_seurat@meta.data$orig.ident)

new_merged_seurat <- merge(
  x = Celltype_Rectum,
  y = list(Celltype_Blood, Celltype_Marrow,Celltype_spleen,Celltype_Stomach,Celltype_Trachea,Celltype_Lymph),
  add.cell.ids = c("Rectum", "Blood", "Marrow","Spleen", "Stomach","Trachea","Lymph")
)

table(new_merged_seurat@meta.data$orig.ident)


#, Celltype_Muscle,Celltype_Skin, Celltype_Stomach, Celltype_spleen, Celltype_SI, Celltype_Marrow, Celltype_Lymph, Celltype_Blood,
file_name <- paste0(celltype,"_merged.rds")
saveRDS(new_merged_seurat,file_name)

Idents(new_merged_seurat) <- new_merged_seurat@meta.data$orig.ident
plot1 <- VlnPlot(new_merged_seurat,features = c("CD19","MS4A1","CD79A"))
plot1



