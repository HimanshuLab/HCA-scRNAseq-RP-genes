library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(EnhancedVolcano)
library(harmony)
library(plotly)


Tissue_name <- "Skin"

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
seurat <- subset(seurat, subset = nFeature_RNA >500 &nCount_RNA >1000  & percent.mt < 25)
seurat

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

# spliting the RNA assay based on patient - information present in metadata
seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$id)
seurat

options(future.globals.maxSize = 1000 * 1024^2)  # Increase to 1000 MiB (1 GB)
#SCTransform after split data normalizes each patient data separately for the genes based on sequencing depths
seurat <- SCTransform(seurat, verbose = TRUE)

#vars.to.regress = "percent.mt" was done for tissues with any negative or high positive corrolation of percent.mt with nFeature_RNA
#seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = TRUE)


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


#nps chosen will be between 20-50. if the elbow shows less that 20 pcs,  50 will be chosen based on variance explained. 
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
p6 <- DimPlot(seurat, reduction = "umap.unintegrated") + #, group.by=c("id")) +
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


#harmony or CCA for batch correction based on patients
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

#Further clustering after batch correction
pcs <- 50
seurat <- FindNeighbors(seurat, reduction = "integrated.cca", dims = 1:pcs)
seurat <- FindClusters(seurat, resolution = 1.2, cluster.name = "CCA_clusters")
seurat <- RunUMAP(seurat, reduction = "integrated.cca", dims = 1:pcs, reduction.name = "umap", metric = "euclidean",seed.use = 10 )
seurat


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

marker <- "marker_genes" #choose marker genes from TableS4

p <- FeaturePlot(seurat, features = marker, reduction = "umap", combine = TRUE)
# Replace 'RNA_snn_res.1' with your specific cluster column name in metadata
p$data$idents <- Idents(seurat)
# Label clusters on the plot
plot <- LabelClusters(plot = p, id = "idents")
plot

vln <- VlnPlot(seurat,features = marker)
vln

p7

file_path <- paste(Tissue_name, "markers",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(plot)
# Close the PNG device
dev.off()

file_path <- paste(Tissue_name, "markers",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(vln)
# Close the PNG device
dev.off()

p7

seurat <- PrepSCTFindMarkers(seurat)
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Select the top 10 markers for each cluster
top50 <- seurat.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)
write.csv(top50,"marrow_222_cluster.csv")

markers <- c("CD3D","TSR1","CA1","TFRC","HBG1","NKG7","GZMK","CD3D","CPA3","KIT","HBA1","S100A8","MPO","HBB","TFRC","SPINK2","PLEK","HBB","CDK2","HBA1","CD3D","HBB","C1QB","ELENA","TPSAB1","MZB1","SPINK2","CD79A","ALAD","IGLL2","JCHAIN","THBS1","HBA1","CD86","CD83","IGFS6","CD1C")

new_idents <- c("Enterocytes","Enterocytes","Enterocytes","CD8T","Enterocytes","Enterocytes","Tuft","Absorptive","Tgammadelta","B","Stem","Macrophage","Plasma","Enteroendocrine","Goblet","Enterocytes","T")
names(new_idents) <- levels(seurat)
seurat <- RenameIdents(seurat, new_idents)
p10 <- DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.2)
p10

seurat@meta.data$celltypes <- Idents(seurat)

#save the processed seurat object
saveRDS(seurat,"SmallIntestine_11_06.rds")
