library(monocle3)
library(SeuratWrappers)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(randomcoloR)
library(ggplot2)
library(ggplot2)
library(dplyr)

setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS")
seurat <- readRDS("Testis_Adult_05_04.rds")

#subsetting seurat for spermatogenesis cells for pseudotime
levels(seurat)

#subsetting the seurat just for 
spermatogenesis <- subset(seurat, idents = c("SSC","DS","Early PS","Late PS","Round Sptd","Elongated Sptd","Sperm"))
spermatogenesis

p10 <- DimPlot(spermatogenesis, reduction = "umap")
ggplotly(p10)

#monocle3 create cell data object based on seurat and since the batch correction and dimensionality reduction is already done it need not be done again
cds <- as.cell_data_set(spermatogenesis, assay = "SCT")
cds

#lower resolution to prevent over clustering and keep a single partition
cds <- cluster_cells(cds, resolution = 0.001)
cds@colData$clusters <- clusters(cds)
p1 <- plot_cells(cds, color_cells_by = "Patient_Replicate", 
                 show_trajectory_graph = FALSE, label_cell_groups = FALSE)
p1


cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)

#this function helps in identifying the cell type that appears the rliest
get_earliest_principal_node <- function(cds, cell_type = "SSC") {
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
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE)
p1 <- plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = FALSE)
p1
saveRDS(cds,"Testis_pseudotime.rds")