library(ComplexHeatmap)
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)

seurat <- readRDS("Plasma_merged.rds")
Idents(seurat) <- seurat@meta.data$orig.ident
seurat <- SCTransform(seurat)
seurat <- PrepSCTFindMarkers(seurat)
Tissue_name <- "Plasma"

#nps chosen will be between 30-50. if the elbow shows less that 30 pcs, 30 will be chosen. 
seurat <- RunPCA(seurat,npcs = 100, verbose = TRUE)
total_variance <- seurat@reductions$pca@misc$total.variance

#calculate the ratio of stdev for each pc with the sum of all variances to understand the percentage variance explained by each PC
std_pca <- seurat[["pca"]]@stdev
var_pca <- std_pca^2
pct <- (var_pca) / total_variance * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
cumu

# Determine the difference between variation of PC and subsequent PC
elbow_pc <- sort(which(((pct[1:length(pct) - 1])^2 - (pct[2:length(pct)])^2) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
elbow_pc

p5 <- ElbowPlot(seurat, ndims = 100)
p5

# choosing the number of PCs
pcs <- 50

# seurat find neighbours to construct SNN graph to apply louvian algorithm on - k.param - Defines k for the k-nearest neighbor algorithm
#As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type.
# k.param value depends on the dataset size and cell types. For datasets with a large number of cells or significant heterogeneity, a higher k.param (e.g., 20-30) may be beneficial to capture broader relationships among cells. Conversely, for smaller datasets or more homogeneous populations, a lower value (e.g., 5-15) might suffice
# kparam values range between 10-50
seurat <- FindNeighbors(seurat, dims = 1:pcs, reduction = "pca")

#Adjust the resolution parameter based on number of clusters
seurat <- FindClusters(seurat, resolution = 1, cluster.name = "unintegrated_clusters")

#plotting UMAP based on louvian clusters using same n_neighbour as k.param in find neighbours. metric can be euclidean or corrolation
seurat <- RunUMAP(seurat, dims = 1:pcs, reduction = "pca", reduction.name = "umap.unintegrated",metric = "correlation",seed.use = 10)

# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
p6 <- DimPlot(seurat, reduction = "umap.unintegrated", group.by=c("orig.ident")) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 22)) 
p6

Idents(seurat) <- seurat@meta.data$orig.ident

#finding rp markers across tissues with tissues as idents
markers <- FindAllMarkers(
  seurat,
  assay = "SCT",
  slot = "data",
  verbose = TRUE
)

rpgenes <- c("RPL39","RPL10","RPS2","RPL21","RPL18A","RPL13A","RPS17","RPS19","RPS18","RPS14","RPL34","RPS12","RPS27A","RPLP1","RPS8","RPL12","RPL23A","RPL18","RPS28","RPSA","RPS15A","RPS3","RPS21","RPL10A","RPL5","RPL24","RPL6","RPL9","RPL22","RPL14","RPL37","RPS16","RPL35A","RPL11","RPS7","RPS23","RPL26","RPL41","RPL7","RPL23","RPS20","RPL27","RPL31","RPL27A","RPS27","RPL3","RPLP0","RPS3A","RPS5","RPL4","RPS29","RPS11","RPL38","RPL37A","RPS9","RPL7A","RPL8","RPS6","RPL32","RPL13","RPS13","RPL36","RPLP2","RPS24","RPL15","RPL19","RPS15","RPL35","RPL30","RPS25","RPL28","RPS10","RPL29","RPS26","RPL17","RPL36A")

plot <- EnhancedVolcano(markers, 
                        rownames(markers),
                        x ="avg_log2FC",
                        y ="p_val_adj",
                        selectLab = rpgenes,
                        subtitle = 'Highlighting RP genes')


plot

filtered_markers <- markers[markers$p_val_adj < 0.05 & (markers$avg_log2FC > 1 | markers$avg_log2FC < -1) ,]

# Filter the markers based on the gene list
subset_markers <- filtered_markers[filtered_markers$gene %in% rpgenes, ]

# Check if any genes exist in the markers
if (nrow(subset_markers) == 0) {
  message("No RP genes from the list exist in the markers.")
} else {
  # Save the subsetted markers to a CSV file
  filename <- paste0(Tissue_name,"RP_DEGs.csv")
  write.csv(subset_markers, filename)
  message(paste("Saved subsetted markers to:", filename))
}

Celltype <- "Stomach_plasma"
markers <- FindMarkers(
  seurat,
  ident.1 = Celltype,   # Change to your cluster/ident
  assay = "SCT",
  slot = "data",
  verbose = TRUE
)

markers$gene <- rownames(markers)

#specific tissue volcano plot
Tissue_name <- "Plasma"
name <- paste0("Stomach"," Vs other tissues")

plot <- EnhancedVolcano(
  markers,
  lab = markers$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  selectLab = rpgenes,
  title = Tissue_name,
  subtitle = name# only label RP genes present in DEGs
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




