library(Seurat)

#infant testis
#when input data is text file
my_data <- read.delim("GSE120506_infant_combined_UMI.txt",skipNul=T)
rownames(my_data) <- my_data[,1]
my_data <- my_data[,-1]

# metadata input from the paper
metadata_file <- read.csv("testis_infant_metadata.csv")
#creating a seurat object using the data input
seurat <- CreateSeuratObject(counts = my_data,assay="RNA", min.cells = 3)
seurat

#add the metadata columns to the seurat object
seurat <- AddMetaData(seurat, metadata = metadata_file, col.name = c("Donor","Replicate","Donor_Replicate","ID"))
# extracting the seurat metadata after adding the patient details to check for the correctness
metadata_seurat <- seurat@meta.data
# to check if the metadata has been added correctly
table(metadata_seurat$ID == rownames(metadata_seurat))
#check if the order of metadata added is right
sum(rownames(metadata_seurat) == metadata_seurat$ID) == nrow(metadata_seurat)


#Adult testis
#when input data is text file
my_data_n <- read.delim("GSE112013_Combined_UMI_table.txt",skipNul=T)
rownames(my_data_n) <- my_data_n[,1]
my_data_n <- my_data_n[,-1]

# metadata input from the paper
metadata_file_n <- read.csv("metadata_testis.csv")

seurat <- CreateSeuratObject(counts = my_data_n,assay="RNA", min.cells = 3)
seurat

#add the metadata columns to the seurat object
seurat <- AddMetaData(seurat, metadata = metadata_file_n, col.name = c("Donor","Replicate","Donor_Replicate","ID"))
# extracting the seurat metadata after adding the petient details to check for the correctness
metadata_seurat <- seurat@meta.data
# to check if the metadata has been added correctly
table(metadata_seurat$ID == rownames(metadata_seurat))


#GSE222368 
sceasy::convertFormat("heme_velo_aggr_5exp_211012.h5ad", from="anndata", to="seurat",
                      outFile='GSE222368.rds')
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


#GSE159929 data input
#blood
data <-read.csv(file = "GSM4850587_Skin_Counts.csv", header = FALSE, sep = ",")
colnames(data) <- data[1,]
rownames(data) <- data[,1]
data <- data[-1,-1]
seurat <- CreateSeuratObject(counts = data, project = "skin", min.cells = 3)

#GSE169062 infertility and adult testis merge

#infertile testis
#when input data is text file
my_data <- read.delim("GSE169062_GeneExpressionMatrix_Infertile.txt",skipNul=T)
rownames(my_data) <- my_data[,1]
my_data <- my_data[,-1]

# metadata input from the paper
metadata_file <- read.csv("GSE169062_metadata.csv")
#creating a seurat object using the data input
seurat_inf <- CreateSeuratObject(counts = my_data,assay="RNA", min.cells = 3)
seurat_inf

#add the metadata columns to the seurat object
seurat_inf <- AddMetaData(seurat_inf, metadata = metadata_file, col.name = c("Donor","Replicate","Donor_Replicate","ID"))
# extracting the seurat metadata after adding the patient details to check for the correctness
metadata_seurat <- seurat_inf@meta.data
# to check if the metadata has been added correctly
table(metadata_seurat$ID == rownames(metadata_seurat))
#check if the order of metadata added is right
sum(rownames(metadata_seurat) == metadata_seurat$ID) == nrow(metadata_seurat)



#Adult testis
#when input data is text file
my_data_n <- read.delim("GSE112013_Combined_UMI_table.txt",skipNul=T)
rownames(my_data_n) <- my_data_n[,1]
my_data_n <- my_data_n[,-1]

# metadata input from the paper
metadata_file_n <- read.csv("metadata_testis.csv")

seurat_n <- CreateSeuratObject(counts = my_data_n,assay="RNA", min.cells = 3)
seurat_n

#add the metadata columns to the seurat object
seurat_n <- AddMetaData(seurat_n, metadata = metadata_file_n, col.name = c("Donor","Replicate","Donor_Replicate","ID"))
# extracting the seurat metadata after adding the petient details to check for the correctness
metadata_seurat <- seurat_n@meta.data
# to check if the metadata has been added correctly
table(metadata_seurat$ID == rownames(metadata_seurat))

#check if the order of metadata added is right
sum(rownames(metadata_seurat) == metadata_seurat$ID) == nrow(metadata_seurat)
seurat <- merge(seurat_inf,seurat_n)
