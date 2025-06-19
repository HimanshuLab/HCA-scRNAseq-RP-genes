library(readr) 
library(Seurat)
library(monocle3)
library(dplyr)


setwd("/Users/aishwaryasharan/Desktop/ScRNA_newanalysis/FinalRDS/RDS_files")

seurat <- readRDS("Testis_Adult_05_04.rds")

#subset to spermatogenesis celltypes
spermatogenesis <- subset(seurat, subset = cluster %in% c("SSC","DS","Early PS","Late PS","Round Sptd","Elongated Sptd","Sperm"))
rpgenes <- c("RPL23","RPS17","RPL36A","RPS20","RPL27A","RPL31","RPL7","RPL27","RPL4","RPS11","RPL38","RPS29","RPL13A","RPSA","RPL10A","RPL23A","RPS2","RPL21","RPS21","RPL18A","RPL18","RPL12","RPS8","RPS28","RPL41","RPL24","RPL9","RPS16","RPL22","RPL17","RPL14","RPS18","RPL6","RPS3","RPS9","RPL35","RPL3","RPS6","RPL5","RPS5","RPLP0","RPS10","RPS19","RPS25","RPS3A","RPS7","RPL26","RPS14","RPS23","RPL37","RPL36","RPL35A","RPS13","RPS27A","RPS12","RPL11","RPL34","RPS26","RPS15A","RPLP1","RPL37A","RPL15","RPL29","RPL7A","RPLP2","RPL8","RPS24","RPS27","RPL30","RPL32","RPS15","RPL19","RPL13","RPL28","RPL39","RPL10")

#extract the patient specific data
patient1 <- subset(spermatogenesis, subset = Patient_Replicate %in% c("Donor1_1", "Donor1_2"))
patient2 <- subset(spermatogenesis, subset = Patient_Replicate %in% c("Donor2_1", "Donor2_2"))
patient3 <- subset(spermatogenesis, subset = Patient_Replicate %in% c("Donor3_1", "Donor3_2"))

#get gene expression data
expr_data <- FetchData(spermatogenesis, vars = rpgenes, assay = "SCT", slot = "data")

# Add identity information
expr_data$ident <- Idents(spermatogenesis)

# Compute average expression per ident
avg_expr <- aggregate(. ~ ident, data = expr_data, FUN = mean)

# Convert to a matrix (genes as rows, idents as columns)
avg_expr_matrix <- as.matrix(avg_expr[,-1])  # Remove the ident column
rownames(avg_expr_matrix) <- avg_expr$ident  # Set row names as idents

# Transpose to get genes as rows and idents as columns
avg_expr_matrix <- t(avg_expr_matrix)

# Print matrix
print(avg_expr_matrix)

#order the celltypes in the right order of spermatogenesis
custom_order <- c("SSC","DS","Early PS","Late PS","Round Sptd","Elongated Sptd","Sperm")  # Define desired order
ordered_matrix <- avg_expr_matrix[, custom_order]
ordered_matrix <- ordered_matrix[rowSums(ordered_matrix) > 0, ]
write.csv(ordered_matrix, "testis_Overall_rp.csv")

# Rename columns if they exactly match original_names
original_names <- c("SSC","DS","Early PS","Late PS","Round Sptd","Elongated Sptd","Sperm")
if (all(colnames(ordered_matrix) == original_names)) {
  print(paste("The colnames are right"))
  colnames(ordered_matrix) <- as.character(1:7)
}

Tissue_name <- "Testis_adult_Overall"
# Convert to data frame and save as TSV
df <- as.data.frame(ordered_matrix)
df <- cbind(RowNames = rownames(ordered_matrix), df)
filename <- paste0(Tissue_name, "_rp.tsv")
write_tsv(df, file = filename)

