library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(EnhancedVolcano)
library(harmony)
library(plotly)


seurat <- readRDS("GSE222368.rds")
# Calculate nFeature_RNA: number of genes with non-zero counts per cell
seurat$nFeature_RNA <- Matrix::colSums(seurat@assays$RNA@counts > 0)

# Calculate nCount_RNA: total counts (UMIs) per cell
seurat$nCount_RNA <- Matrix::colSums(seurat@assays$RNA@counts)

# Calculate the number of cells each gene is detected in (non-zero counts)
gene_cell_counts <- Matrix::rowSums(seurat@assays$RNA@counts > 0)

# Identify genes expressed in at least 20 cells
genes_to_keep <- names(gene_cell_counts[gene_cell_counts >= 20])

# Subset the Seurat object to keep only these genes
seurat_filtered <- subset(seurat, features = genes_to_keep)

# Optional: Check the number of genes before and after filtering
cat("Number of genes before filtering:", nrow(seurat@assays$RNA@counts), "\n")
cat("Number of genes after filtering:", nrow(seurat_filtered@assays$RNA@counts), "\n")


Tissue_name <- "Marrow_222"

#The input for this file is from the output of data_input file
seurat <- seurat_filtered
# Identify the percentage of mt genes in each cells and add it to metadata
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
# these tell us about the distribution of the data
p1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1


file_path <- paste(Tissue_name, "percent.mt_Vlnplot_before",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(p1)
# Close the PNG device
dev.off()


seurat

#subsetting the data based on the cutoffs mentioned in the respective studies, thus change the parameters accordingly
#testis
#nFeatureRNA- total number of genes expressed in the cell, nCountsRNA - total numbe rof molecules expressed in the cell
seurat <- subset(seurat, subset = nFeature_RNA >200  & percent.mt < 15)
seurat1 <- seurat

FeatureScatter(seurat, feature1 = "percent.mt", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2


file_path <- paste(Tissue_name, "percent.mt_Vlnplot_after",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(p2)
# Close the PNG device
dev.off()

#seurat <- JoinLayers(seurat) #use this if this is a merged seurat object

seurat <- subset(seurat, subset = batch %in% c("D1d0a","D1d0b","D1d3b","D1d6a","D1d6b","D2d0a","D2d0b","D2d3b","D2d6a","D2d6b","N1d0a","N1d6a","N2d0a","N2d0b","N2d3b","N2d6a","N2d6b","N3d0b","N3d3b","N3d6b"))
# spliting the RNA assay based on patient - information present in metadata
seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$id)
seurat

FeatureScatter(seurat, feature1 = "percent.mt", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

options(future.globals.maxSize = 1000 * 1024^2)  # Increase to 1000 MiB (1 GB)
#SCTransform after split data normalizes each patient data seperatly for the genes based on sequencing depths
#seurat <- SCTransform(seurat, verbose = TRUE)
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = TRUE)


#checking if the sctransformed count has similar depth for all cells
p3 <- hist(seurat@meta.data$nCount_SCT , breaks = 100)
p3

file_path <- paste(Tissue_name, "hist_afterSCT",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p3)
# Close the PNG device
dev.off()

p4 <- hist(seurat@meta.data$nCount_RNA , breaks = 100)
p4
file_path <- paste(Tissue_name, "hist_beforeSCT",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p4)
# Close the PNG device
dev.off()

#checking if the median is similar for both raw and corrected counts
median(seurat@meta.data$nCount_SCT)
median(seurat@meta.data$nCount_RNA)


#nps chosen will be between 20-50. if the elbow shows less that 20 pcs, 50 will be chosen based on variance explained. 
seurat <- RunPCA(seurat,npcs = 100, verbose = TRUE)
total_variance <- seurat@reductions$pca@misc$total.variance

#calculate the ratio of stdev for each pc with the sum of all variances to understand the percentage variance explained by each PC
std_pca <- seurat[["pca"]]@stdev
var_pca <- std_pca^2
pct <- (var_pca) / total_variance * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
cumu

file_name <- paste(Tissue_name, "PCA_varexplained",".csv", sep = "")
write.csv(cumu,file_name)

# Determine the difference between variation of PC and subsequent PC
elbow_pc <- sort(which(((pct[1:length(pct) - 1])^2 - (pct[2:length(pct)])^2) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
elbow_pc

p5 <- ElbowPlot(seurat, ndims = 100)
p5
file_path <- paste(Tissue_name, "elbowplot",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p5)
# Close the PNG device
dev.off()


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
seurat <- RunUMAP(seurat, dims = 1:pcs, reduction = "pca", reduction.name = "umap.unintegrated",metric = "euclidean",seed.use = 10)

# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
p6 <- DimPlot(seurat, reduction = "umap.unintegrated", group.by=c("id")) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 22)) + # Increase y-axis label font size
  labs(title = "Before Integration") 
p6

file_path <- paste(Tissue_name, "umap_preinteg",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p6)
# Close the PNG device
dev.off()


#harmony for batch correction based on patients
seurat <- RunHarmony(
  seurat,
  group.by.vars = "id",  # Specify the grouping variable
  reduction.use = "pca",# Use PCA as input for batch correction
  plot_convergence = TRUE,
  dims.use = 1:100 #generally all Pcs are used
)

seurat <- IntegrateLayers(seurat, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "SCT",
  verbose = TRUE
)

pcs <- 50
seurat <- FindNeighbors(seurat, reduction = "integrated.cca", dims = 1:pcs)
seurat <- FindClusters(seurat, resolution = 1.2, cluster.name = "CCA_clusters")
seurat <- RunUMAP(seurat, reduction = "integrated.cca", dims = 1:pcs, reduction.name = "umap", metric = "euclidean",seed.use = 10 )
#seurat <- RunTSNE(seurat, reduction = "harmony", dims = 1:pcs, reduction.name = "tsne", metric = "correlation",seed.use = 345, verbose = TRUE)
seurat

seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:pcs)
seurat <- FindClusters(seurat, resolution = 1.2, cluster.name = "harmony_clusters")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:pcs, reduction.name = "umap_harmony", metric = "correlation",seed.use = 10 )

# #checking different seed s for a different distribution
p7 <- DimPlot(seurat, reduction = "umap", group.by = c("CCA_clusters"), label = TRUE) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 20)) + # Increase y-axis label font size
  labs(title = "After Integration") 


p7

p8 <- DimPlot(seurat, reduction = "umap", group.by = c("id"),pt.size = 0.01) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 22)) + # Increase y-axis label font size
  labs(title = "After Integration") 

p7 + p8

file_path <- paste(Tissue_name, "umap_postinteg_clusters",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p7)
# Close the PNG device
dev.off()

file_path <- paste(Tissue_name, "umap_postinteg",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p8)
# Close the PNG device
dev.off()

file_name <- paste(Tissue_name, "umap_pre-post_integ",".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(p6+p8)
dev.off()

p9 <- DoHeatmap(seurat, features = allmarkers)
p9

file_path <- paste(Tissue_name, "cluster_marker_heatmap",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(p9)
# Close the PNG device
dev.off()


p <- FeaturePlot(seurat, features = Fibro, reduction = "umap", combine = TRUE)
# Replace 'RNA_snn_res.1' with your specific cluster column name in metadata
p$data$idents <- Idents(seurat)
# Label clusters on the plot
plot <- LabelClusters(plot = p, id = "idents")
plot

vln <- VlnPlot(seurat,features = "KIT")
vln

p7

file_path <- paste(Tissue_name, "Fibro",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(plot)
# Close the PNG device
dev.off()

file_path <- paste(Tissue_name, "Fibro",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(vln)
# Close the PNG device
dev.off()

p7

seurat1 <- seurat
seurat1 <- PrepSCTFindMarkers(seurat1)
seurat.markers <- FindAllMarkers(seurat1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Select the top 10 markers for each cluster
top50 <- seurat.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)
write.csv(top50,"marrow_222_cluster.csv")

markers <- c("CD3D","TSR1","CA1","TFRC","HBG1","NKG7","GZMK","CD3D","CPA3","KIT","HBA1","S100A8","MPO","HBB","TFRC","SPINK2","PLEK","HBB","CDK2","HBA1","CD3D","HBB","C1QB","ELENA","TPSAB1","MZB1","SPINK2","CD79A","ALAD","IGLL2","JCHAIN","THBS1","HBA1","CD86","CD83","IGFS6","CD1C")
VlnPlot(seurat,features = markers)
new_idents <- c("Ery_pro","Erythrocyte1","Erythrocyte1","Erythrocyte1","Tcells","Erythrocyte1","Tcells","Erythrocyte2","Erythrocyte2","Erythrocyte2","Erythrocyte1","Erythrocyte1","CLP","Tcells","NKT CD8+","Erythrocyte1","Neutrophil","Erythrocyte2","Macrophage/mono","Monocytes","Erythrocyte1","unkown","FLT3+ CLP","Plasmocytoid dendritic","Erythrocyte1","Erythrocyte2","DC FLT3+","Erythrocyte2","B","Erythrocyte1","Erythrocyte2","Megakaryocyte","Erythrocyte2","DC")
names(new_idents) <- levels(seurat)
seurat <- RenameIdents(seurat, new_idents)
p10 <- DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.2)
p10


seurat@meta.data$celltypes <- Idents(seurat)
#seurat<- subset(seurat,subset = celltypes != "unknown" )
saveRDS(seurat,"marrow_14_6.rds")
