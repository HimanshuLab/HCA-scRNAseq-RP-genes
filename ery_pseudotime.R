library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(randomcoloR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)

seurat <- readRDS("marrow_14_6.rds")
levels(seurat)

ery <- subset(seurat, idents = c("Erythrocyte1","Erythrocyte2","Ery_pro"))
saveRDS(ery,"GSE222_ery_unclustered.rds")

#confirm erythro and progenitor markers in this 
#----
VlnPlot(ery, features = c("KIT","CA2","TFRC","HBB","HBA1","CD36","GATA1"))

Tissue_name <- "Erythrocytes_222"
erythro <- readRDS("GSE222_ery_unclustered.rds")
table(ery@meta.data$predicted_doublet)
erythro <- subset(erythro, subset = predicted_doublet %in% "FALSE")
table(erythro@meta.data$batch)

erythro <- JoinLayers(erythro, assay = "RNA")

#creating a new seurat object and doing UMAP again just for erythrocytes
ery <- CreateSeuratObject(counts = erythro@assays$RNA$counts)
cols_to_copy <- c("day", "barcode","id","batch","doublet_score","predicted_doublet","prob","trajectoryb","hap","exp")
ery@meta.data[, cols_to_copy] <- erythro@meta.data[, cols_to_copy]
table(ery@meta.data$exp)
ery[["percent.mt"]] <- PercentageFeatureSet(ery, pattern = "^MT-")
ery[["RNA"]] <- split(ery[["RNA"]], f = ery$id)
ery

FeatureScatter(ery, feature1 = "percent.mt", feature2 = "nFeature_RNA")
FeatureScatter(ery, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

options(future.globals.maxSize = 8 * 1024^3)  # Increase to 1000 MiB (1 GB)
#SCTransform after split data normalizes each patient data seperatly for the genes based on sequencing depths
ery <- SCTransform(ery, verbose = TRUE)


Sys.time()
#nps chosen will be between 20-50. if the elbow shows less that 20 pcs, 20 will be chosen. 
ery <- RunPCA(ery,npcs = 100, verbose = TRUE)
Sys.time()


total_variance <- ery@reductions$pca@misc$total.variance

#calculate the ratio of stdev for each pc with the sum of all variances to understand the percentage variance explained by each PC
std_pca <- ery[["pca"]]@stdev
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

Tissue_name <- "erythro_222"
p5 <- ElbowPlot(ery, ndims = 100)
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

# ery find neighbours to construct SNN graph to apply louvian algorithm on - k.param - Defines k for the k-nearest neighbor algorithm
#As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type.
# k.param value depends on the dataset size and cell types. For datasets with a large number of cells or significant heterogeneity, a higher k.param (e.g., 20-30) may be beneficial to capture broader relationships among cells. Conversely, for smaller datasets or more homogeneous populations, a lower value (e.g., 5-15) might suffice
# kparam values range between 10-50
ery <- FindNeighbors(ery, dims = 1:pcs, reduction = "pca")

#Adjust the resolution parameter based on number of clusters
ery <- FindClusters(ery, resolution = 1, cluster.name = "unintegrated_clusters")

#plotting UMAP based on louvian clusters using same n_neighbour as k.param in find neighbours. metric can be euclidean or corrolation
ery <- RunUMAP(ery, dims = 1:pcs, reduction = "pca", reduction.name = "umap.unintegrated",metric = "euclidean",seed.use = 10)

# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
p6 <- DimPlot(ery, reduction = "umap.unintegrated", group.by=c("id")) +
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


library(harmony)
#harmony for batch correction based on patients
ery <- RunHarmony(
  ery,
  group.by.vars = "id",  # Specify the grouping variable
  reduction.use = "pca",# Use PCA as input for batch correction
  plot_convergence = TRUE,
  dims.use = 1:100 #generally all Pcs are used
)

pcs <- 50
ery <- FindNeighbors(ery, reduction = "harmony", dims = 1:pcs)
ery <- FindClusters(ery, resolution = 1.2, cluster.name = "harmony_clusters", algorithm = 1, verbose = TRUE)
ery <- RunUMAP(ery, reduction = "harmony", dims = 1:pcs, reduction.name = "umap", metric = "euclidean",seed.use = 10)

# #checking different seed s for a different distribution
p7 <- DimPlot(seurat, reduction = "umap", label = TRUE) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 20)) + # Increase y-axis label font size
  labs(title = "After Integration") 
p7


p8 <- DimPlot(ery, reduction = "umap", group.by = c("id"),pt.size = 0.01) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 22)) + # Increase y-axis label font size
  labs(title = "After Integration") 
p8

VlnPlot(ery, features = c("KIT"))

p <- FeaturePlot(ery, features = "KIT", reduction = "umap", combine = TRUE)
# Replace 'RNA_snn_res.1' with your specific cluster column name in metadata
p$data$idents <- Idents(ery)
# Label clusters on the plot
plot <- LabelClusters(plot = p, id = "idents")
plot
#----

#finding differentially expressed markers in ery
#----
ery1 <- ery
ery1 <- PrepSCTFindMarkers(ery1)
ery.markers <- FindAllMarkers(ery1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top50 <- ery.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)
write.csv(ery.markers,"ery_222_cluster.csv")
#----

# Assigning cell type
#----
#based on findallmarkers and Violin plot, cluster 1,17 has highest kit expression
new_idents <- c("0","Progenitor","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","Progenitor","18","19","20","21","22","23","24","25")
names(new_idents) <- levels(ery)
ery <- RenameIdents(ery, new_idents)
p10 <- DimPlot(ery, reduction = "umap", label = TRUE, pt.size = 0.2)
p10
saveRDS(ery,"erythrocytes_222_n1n2n3d1d2_12_6_noreg.rds")
#----

#pseudotime
#----
ery <- readRDS("erythrocytes_222_n1n2n3d1d2_12_6_noreg.rds")
cds <- as.cell_data_set(ery, assay = "SCT")
cds

levels(ery)
#lower resolution to prevent over clustering and keep a single partition
cds <- cluster_cells(cds, resolution = 0.0001)
cds@colData$clusters <- clusters(cds)
p1 <- plot_cells(cds, color_cells_by = "clusters", 
                 show_trajectory_graph = FALSE, label_cell_groups = FALSE)
p1

# Create a vector of 1s for all cells
single_partition <- rep(1, ncol(cds))

# Assign cell names to the partition vector
names(single_partition) <- colnames(cds)

# Convert to factor
single_partition <- as.factor(single_partition)

# Assign to the partition slot in the CDS object
cds@clusters$UMAP$partitions <- single_partition

plot_cells(cds, color_cells_by = 'partition', label_cell_groups = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)


get_earliest_principal_node <- function(cds, cell_type = "Progenitor") {
  cell_ids <- which(colData(cds)[, "ident"] == cell_type)
  
  closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)$UMAP)$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,
  ]))))]
  
  root_pr_nodes
}

get_earliest_principal_node(cds)

#ordering cells based on pseudotime
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
plot_cells(cds_D, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = TRUE)

#cds <- readRDS("ery_cds_pseudotime_n1n2n3d1d2.rds")
saveRDS(cds,"ery_cds_pseudotime_n1n2n3d1d2.rds")


# # Extract pseudotime values from the Monocle3 object
pseudotime_values <- pseudotime(cds)
# Ensure the cell names in `cds_N` match those in the Seurat object
cell_names <- colnames(cds)

if (identical(cell_names, colnames(ery))) {
  # Add pseudotime to Seurat metadata
  ery@meta.data$pseudotime <- pseudotime_values[cell_names]
} else {
  stop("Cell names do not match between cds and Seurat object.")
}

saveRDS(ery,"erythrocytes_222_n1n2n3d1d2_16_6_assigned.rds")
#----

#pseudotime plots
#----
ery@meta.data$progenitor <- Idents(ery)
Idents(ery) <- ery@meta.data$stage

# Subset for N1-N3 group
cds_N <- cds[, colData(cds)$id %in% paste0("N", 1:3)]
# Subset for D1 group
cds_D1 <- cds[, colData(cds)$id %in% c("D1")]
# Subset for D2 group
cds_D2 <- cds[, colData(cds)$id %in% c("D2")]


plot <- plot_cells(cds_N, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                   label_branch_points = FALSE)
file_name <- paste("pseudotime_normal",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("pseudotime_normal",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()

plot <- plot_cells(cds_D1, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                   label_branch_points = FALSE)

file_name <- paste("pseudotime_DBA1",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("pseudotime_DBA1",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()

plot <- plot_cells(cds_D2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                   label_branch_points = FALSE)

file_name <- paste("pseudotime_DBA2",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("pseudotime_DBA2",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()
#----

#Input for DPGP
#----
ery <- readRDS("erythrocytes_222_n1n2n3d1d2_16_6_assigned.rds")
cds <- readRDS("ery_cds_pseudotime_n1n2n3d1d2_16_6.rds")
pseudotime <- ery@meta.data$pseudotime  # Verify column name[2][4]
expr_mat <- t(expression_matrix)

#assigning gene names as a row
rpgenes <- c("RPL23","RPS17","RPL36A","RPS20","RPL27A","RPL31","RPL7","RPL27","RPL4","RPS11","RPL38","RPS29","RPL13A","RPSA","RPL10A","RPL23A","RPS2","RPL21","RPS21","RPL18A","RPL18","RPL12","RPS8","RPS28","RPL41","RPL24","RPL9","RPS16","RPL22","RPL17","RPL14","RPS18","RPL6","RPS3","RPS9","RPL35","RPL3","RPS6","RPL5","RPS5","RPLP0","RPS10","RPS19","RPS25","RPS3A","RPS7","RPL26","RPS14","RPS23","RPL37","RPL36","RPL35A","RPS13","RPS27A","RPS12","RPL11","RPL34","RPS26","RPS15A","RPLP1","RPL37A","RPL15","RPL29","RPL7A","RPLP2","RPL8","RPS24","RPS27","RPL30","RPL32","RPS15","RPL19","RPL13","RPL28","RPL39","RPL10")
ery_N <- subset(ery, subset = id %in% c("N1","N2","N3"))
ery_D <- subset(ery, subset = id %in% c("D1","D2"))

#assigning pseudotime as a column
cds_N@colData$pseudotime <- pseudotime(cds_N)
expr_mat <- ery_N[["SCT"]]$data[rpgenes, ] %>% as.matrix()
pseudotime <- pseudotime(cds_N)
table(colnames(expr_mat) == colnames(ery_N))

data <- data.frame(pseudotime = pseudotime, t(expr_mat))

# Sort the data frame by pseudotime
data <- data[order(data$pseudotime), ]

# Get the breaks used by cut() for normal and apply it for all
breaks <- seq(min(data$pseudotime, na.rm = TRUE), max(data$pseudotime, na.rm = TRUE), length.out = 8)

# Create bins with these breaks and labels
time_bins <- cut(data$pseudotime, breaks = breaks, labels = paste0("Stage_", 1:7), include.lowest = TRUE)

bin_ranges <- data.frame(
  Stage = paste0("Stage_", 1:7),
  Lower = breaks[-length(breaks)],
  Upper = breaks[-1]
)

print(bin_ranges)

data$time_bin <- time_bins

# Extract sorted expression matrix and pseudotime
sorted_expression_matrix <- data[, -1]  # Remove pseudotime column
sorted_expression_matrix <- sorted_expression_matrix[, -ncol(sorted_expression_matrix)]
colnames(sorted_expression_matrix)
table(data$time_bin)

agg_expression <- aggregate(sorted_expression_matrix,
                            by = list(time_bins),
                            FUN = mean)

# Rename columns for clarity
colnames(agg_expression)[1] <- "Time_Interval"

write.csv(t(agg_expression),"ery_normal_pseudotime_aggregate_7bin.csv")

avg_expr <- as.data.frame(t(agg_expression))
colnames(avg_expr) <- avg_expr[1,]
avg_expr <- avg_expr[-1,]

original_names <- c("Stage_1","Stage_2","Stage_3","Stage_4","Stage_5","Stage_6","Stage_7")
ordered_matrix <- avg_expr[, original_names]

# Rename columns for clarity
colnames(agg_expression)[1] <- "Time_Interval"

write.csv(ordered_matrix,"ery_Normal_pseudotime.csv")

original_names <- c("Stage_1","Stage_2","Stage_3","Stage_4","Stage_5","Stage_6","Stage_7")
if (all(colnames(ordered_matrix) == original_names)) {
  print(paste("The colnames are right"))
  colnames(ordered_matrix) <- as.character(1:7)
}

Tissue_name <- "Erythro_Normal_overall"
# Convert to data frame and save as TSV
df <- as.data.frame(ordered_matrix)
df <- cbind(RowNames = rownames(ordered_matrix), df)
filename <- paste0(Tissue_name, "_rp.tsv")
write_tsv(df, file = filename)

#assigning pseudotime as a column
cds_D@colData$pseudotime <- pseudotime(cds_D)

# # Extract pseudotime values from the Monocle3 object
pseudotime_values <- pseudotime(cds_D)
# Ensure the cell names in `cds_N` match those in the Seurat object
cell_names <- colnames(cds_D)

if (all(cell_names %in% colnames(ery_D))) {
  # Add pseudotime to Seurat metadata
  ery_D@meta.data$pseudotime <- pseudotime_values[cell_names]
} else {
  stop("Cell names do not match between cds and Seurat object.")
}

# for individual donors runt hsi for loop
donor_list <- c("N1","N2","N3")
donor_list <- c("D1","D2")
for (donor_name in donor_list) { 
  donor <- subset(ery_D, subset = id == donor_name)
  expr_data <- FetchData(donor, vars = rpgenes, assay = "SCT", slot = "data")
  pseudotime <- donor@meta.data$pseudotime
  data <- data.frame(pseudotime = pseudotime, expr_data)
  
  # Sort the data frame by pseudotime
  data <- data[order(data$pseudotime), ]
  
  # Create bins with these breaks and labels
  time_bins <- cut(data$pseudotime, breaks = breaks, labels = paste0("Stage_", 1:7), include.lowest = TRUE)

  bin_ranges <- data.frame(
    Stage = paste0("Stage_", 1:7),
    Lower = breaks[-length(breaks)],
    Upper = breaks[-1]
  )
  
  print(bin_ranges)
  
  data$time_bin <- time_bins
  print(donor_name)
  print(table(data$time_bin))
  
  # Extract sorted expression matrix and pseudotime
  sorted_expression_matrix <- data[, -1]  # Remove pseudotime column
  sorted_expression_matrix <- sorted_expression_matrix[, -ncol(sorted_expression_matrix)]
  colnames(sorted_expression_matrix)
  table(data$time_bin)
  
  agg_expression <- aggregate(sorted_expression_matrix,
                              by = list(time_bins),
                              FUN = mean)
  
  avg_expr <- as.data.frame(t(agg_expression))
  colnames(avg_expr) <- avg_expr[1,]
  avg_expr <- avg_expr[-1,]
  
  original_names <- c("Stage_1","Stage_2","Stage_3","Stage_4","Stage_5","Stage_6","Stage_7")
  ordered_matrix <- avg_expr[, original_names]
  
  # Rename columns for clarity
  colnames(agg_expression)[1] <- "Time_Interval"
  
  Tissue_name <- "Erythro_normal_seperate"
  filename <- paste0(Tissue_name,"_",donor_name,".csv")
  write.csv(ordered_matrix,filename)
  
  original_names <- c("Stage_1","Stage_2","Stage_3","Stage_4","Stage_5","Stage_6","Stage_7")
  if (all(colnames(ordered_matrix) == original_names)) {
    print(paste("The colnames are right"))
    colnames(ordered_matrix) <- as.character(1:7)
  }
  
  # Convert to data frame and save as TSV
  df <- as.data.frame(ordered_matrix)
  df <- cbind(RowNames = rownames(ordered_matrix), df)
  filename <- paste0(Tissue_name,"_",donor_name, "_rp.tsv")
  write_tsv(df, file = filename)
}
#----


#Assigning stages to the ery seurat object
#----
pt <- pseudotime(cds)
# Assign bin labels
pt_bins <- cut(pt, breaks = breaks, labels = paste0("Stage_", 1:7), include.lowest = TRUE)
colData(cds)$pt_bin <- pt_bins

# # Extract pseudotime values from the Monocle3 object
stage <- cds@colData$pt_bin
# Ensure the cell names in `cds_N` match those in the Seurat object
cell_names <- colnames(cds)


if (identical(cell_names, colnames(ery))) {
  # Add pseudotime to Seurat metadata
  ery@meta.data$stage <- stage
} else {
  stop("Cell names do not match between cds and Seurat object.")
}

# Subset for N1-N6 group
cds_N <- cds[, colData(cds)$id %in% paste0("N", 1:3)]
# Subset for D1-D2 group
cds_D1 <- cds[, colData(cds)$id %in% c("D1")]
cds_D2 <- cds[, colData(cds)$id %in% c("D2")]


# confirming the correct assignment of stages
p <- plot_cells(cds_N, color_cells_by = "pt_bin",group_label_size = 3)
p <- plot_cells(cds_D1, color_cells_by = "pt_bin",group_label_size = 3)
p <- plot_cells(cds_D2, color_cells_by = "pt_bin",group_label_size = 3)
#----

#marker, and upstream gene expression
#----
ery_N <- subset(ery, subset = id %in% c("N1","N2","N3"))
ery_D1 <- subset(ery, subset = id %in% c("D1"))
ery_D2 <- subset(ery, subset = id %in% c("D2"))

library(ggplot2)
gene_of_interest <- c("KIT","CA2","TFRC","HBB","HBA1","CD36")
#markers_Normal
#----
# Get expression for selected genes (genes x cells)
expression_mat<- ery_N[["SCT"]]$data[gene_of_interest, ] %>% as.matrix()
identical(colnames(expression_mat),colnames(ery_N))
pseudotime <- ery_N@meta.data$pseudotime


# Convert to long format
plot_df <- as.data.frame(as.matrix(expression_mat)) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    -gene,
    names_to = "cell",
    values_to = "expr"
  )

# Create a metadata data.frame with cell names and pt_bin
meta_df <- data.frame(
  cell = colnames(ery_N),
  pt_bin = ery_N@meta.data$stage,
  pseudotime = pseudotime
)

# Join plot_df with meta_df by 'cell'
plot_df <- left_join(plot_df, meta_df, by = "cell")


library(ggplot2)

plot <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = pt_bin )) +
  geom_point(alpha = 2, size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~gene, scales = "fixed") +
  labs(
    x = "Pseudotime",   # Replace with your desired x-axis label
    y = "Normalized expression"    # Replace with your desired y-axis label
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x  = element_text(size = 14),  # X-axis tick label font size
    axis.text.y  = element_text(size = 14)   # Y-axis tick label font size
  )
plot
file_name <- paste("ery_markers_Normal",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("ery_markers_Normal",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()


#----
#markers_DBA1
#----
# Get expression for selected genes (genes x cells)
expression_mat<- ery_D1[["SCT"]]$data[gene_of_interest, ] %>% as.matrix()
identical(colnames(expression_mat),colnames(ery_D1))
pseudotime <- ery_D1@meta.data$pseudotime


# Convert to long format
plot_df <- as.data.frame(as.matrix(expression_mat)) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    -gene,
    names_to = "cell",
    values_to = "expr"
  )

# Create a metadata data.frame with cell names and pt_bin
meta_df <- data.frame(
  cell = colnames(ery_D1),
  pt_bin = ery_D1@meta.data$stage,
  pseudotime = pseudotime
)

# Join plot_df with meta_df by 'cell'
plot_df <- left_join(plot_df, meta_df, by = "cell")


library(ggplot2)

plot <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = pt_bin )) +
  geom_point(alpha = 2, size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~gene, scales = "fixed") +
  labs(
    x = "Pseudotime",   # Replace with your desired x-axis label
    y = "Normalized expression"    # Replace with your desired y-axis label
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x  = element_text(size = 14),  # X-axis tick label font size
    axis.text.y  = element_text(size = 14)   # Y-axis tick label font size
  )
plot
file_name <- paste("ery_markers_DBA1",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("ery_markers_DBA1",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()

#markers_DBA2
#----
# Get expression for selected genes (genes x cells)
expression_mat<- ery_D2[["SCT"]]$data[gene_of_interest, ] %>% as.matrix()
identical(colnames(expression_mat),colnames(ery_D2))
pseudotime <- ery_D2@meta.data$pseudotime


# Convert to long format
plot_df <- as.data.frame(as.matrix(expression_mat)) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    -gene,
    names_to = "cell",
    values_to = "expr"
  )

# Create a metadata data.frame with cell names and pt_bin
meta_df <- data.frame(
  cell = colnames(ery_D2),
  pt_bin = ery_D2@meta.data$stage,
  pseudotime = pseudotime
)

# Join plot_df with meta_df by 'cell'
plot_df <- left_join(plot_df, meta_df, by = "cell")


library(ggplot2)

plot <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = pt_bin )) +
  geom_point(alpha = 2, size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~gene, scales = "fixed") +
  labs(
    x = "Pseudotime",   # Replace with your desired x-axis label
    y = "Normalized expression"    # Replace with your desired y-axis label
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x  = element_text(size = 14),  # X-axis tick label font size
    axis.text.y  = element_text(size = 14)   # Y-axis tick label font size
  )
plot
file_name <- paste("ery_markers_DBA2",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("ery_markers_DBA2",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()

gene_of_interest <- c("GATA1","SPI1","RPS19")
#Upstream regulators_Normal
#----
# Get expression for selected genes (genes x cells)
expression_mat<- ery_N[["SCT"]]$data[gene_of_interest, ] %>% as.matrix()
identical(colnames(expression_mat),colnames(ery_N))
pseudotime <- ery_N@meta.data$pseudotime


# Convert to long format
plot_df <- as.data.frame(as.matrix(expression_mat)) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    -gene,
    names_to = "cell",
    values_to = "expr"
  )

# Create a metadata data.frame with cell names and pt_bin
meta_df <- data.frame(
  cell = colnames(ery_N),
  pt_bin = ery_N@meta.data$stage,
  pseudotime = pseudotime
)

# Join plot_df with meta_df by 'cell'
plot_df <- left_join(plot_df, meta_df, by = "cell")


library(ggplot2)

plot <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = pt_bin )) +
  #geom_point(alpha = 2, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~gene, scales = "fixed") +
  labs(
    x = "Pseudotime",   # Replace with your desired x-axis label
    y = "Normalized expression"    # Replace with your desired y-axis label
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x  = element_text(size = 14),  # X-axis tick label font size
    axis.text.y  = element_text(size = 14)   # Y-axis tick label font size
  )
plot
file_name <- paste("ery_upstream_Normal",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("ery_upstream_Normal",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()

#----
#Upstream regulators_DBA1
#----
# Get expression for selected genes (genes x cells)
expression_mat<- ery_D1[["SCT"]]$data[gene_of_interest, ] %>% as.matrix()
identical(colnames(expression_mat),colnames(ery_D1))
pseudotime <- ery_D1@meta.data$pseudotime


# Convert to long format
plot_df <- as.data.frame(as.matrix(expression_mat)) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    -gene,
    names_to = "cell",
    values_to = "expr"
  )

# Create a metadata data.frame with cell names and pt_bin
meta_df <- data.frame(
  cell = colnames(ery_D1),
  pt_bin = ery_D1@meta.data$stage,
  pseudotime = pseudotime
)

# Join plot_df with meta_df by 'cell'
plot_df <- left_join(plot_df, meta_df, by = "cell")


library(ggplot2)

plot <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = pt_bin )) +
  #geom_point(alpha = 2, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~gene, scales = "fixed") +
  labs(
    x = "Pseudotime",   # Replace with your desired x-axis label
    y = "Normalized expression"    # Replace with your desired y-axis label
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x  = element_text(size = 14),  # X-axis tick label font size
    axis.text.y  = element_text(size = 14)   # Y-axis tick label font size
  )
plot
file_name <- paste("ery_upstream_DBA1",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("ery_upstream_DBA1",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()

#Upstream regulators_DBA2
#----
# Get expression for selected genes (genes x cells)
expression_mat<- ery_D2[["SCT"]]$data[gene_of_interest, ] %>% as.matrix()
identical(colnames(expression_mat),colnames(ery_D2))
pseudotime <- ery_D2@meta.data$pseudotime


# Convert to long format
plot_df <- as.data.frame(as.matrix(expression_mat)) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    -gene,
    names_to = "cell",
    values_to = "expr"
  )

# Create a metadata data.frame with cell names and pt_bin
meta_df <- data.frame(
  cell = colnames(ery_D2),
  pt_bin = ery_D2@meta.data$stage,
  pseudotime = pseudotime
)

# Join plot_df with meta_df by 'cell'
plot_df <- left_join(plot_df, meta_df, by = "cell")


library(ggplot2)

plot <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = pt_bin )) +
  #geom_point(alpha = 2, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~gene, scales = "fixed") +
  labs(
    x = "Pseudotime",   # Replace with your desired x-axis label
    y = "Normalized expression"    # Replace with your desired y-axis label
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x  = element_text(size = 14),  # X-axis tick label font size
    axis.text.y  = element_text(size = 14)   # Y-axis tick label font size
  )
plot
file_name <- paste("ery_upstream_DBA2",".svg", sep = "")
svg(file_name,width = 6, height = 4)
print(plot)
dev.off()

file_path <- paste("ery_upstream_DBA2",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 6, height = 4, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()
#----
