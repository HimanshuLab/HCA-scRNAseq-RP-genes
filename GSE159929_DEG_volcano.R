library(Seurat)
library(EnhancedVolcano)

filename <-"Liver_12_04.rds"
seurat <- readRDS(filename)
table(Idents(seurat))
Tissue_name <- "Liver"

rpgenes <- c("RPL23","RPS17","RPL36A","RPS20","RPL27A","RPL31","RPL7","RPL27","RPL4","RPS11","RPL38","RPS29","RPL13A","RPSA","RPL10A","RPL23A","RPS2","RPL21","RPS21","RPL18A","RPL18","RPL12","RPS8","RPS28","RPL41","RPL24","RPL9","RPS16","RPL22","RPL17","RPL14","RPS18","RPL6","RPS3","RPS9","RPL35","RPL3","RPS6","RPL5","RPS5","RPLP0","RPS10","RPS19","RPS25","RPS3A","RPS7","RPL26","RPS14","RPS23","RPL37","RPL36","RPL35A","RPS13","RPS27A","RPS12","RPL11","RPL34","RPS26","RPS15A","RPLP1","RPL37A","RPL15","RPL29","RPL7A","RPLP2","RPL8","RPS24","RPS27","RPL30","RPL32","RPS15","RPL19","RPL13","RPL28","RPL39","RPL10")
levels(seurat)

Celltype <- "Macrophage"

markers <- FindMarkers(
  seurat,
  ident.1 = Celltype,   # Change to your cluster/ident
  assay = "SCT",
  slot = "data",
  verbose = TRUE
)
markers$gene <- rownames(markers)


name <- paste0(Tissue_name,"- ",Celltype," Vs other celltypes")

plot <- EnhancedVolcano(
  markers,
  lab = markers$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  selectLab = rpgenes,
  title = name # only label RP genes present in DEGs
  )

plot

name1 <- paste0(Tissue_name,"- ",Celltype,"DEG")
file_name <- paste(name1,"heatmap",".svg", sep = "")
svg(file_name,width = 10, height = 8)
plot
dev.off()

file_path <- paste(name1,"heatmap",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot
# Close the PNG device
dev.off()

