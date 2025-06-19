seurat <- readRDS("testis_cntrl_inf_KS_07_04_recent.rds")
levels(seurat)
metadata <- seurat@meta.data

#subsetting KS-2 due to vaery low cell numbers compared to others
seurat <- subset(seurat, subset = DonorCategory != "Klinefelter syndrome-2")

#subset only secondary infertility for germ cells
# seurat <- subset(seurat, subset = DonorCategory != "Klinefelter syndrome-1")
# seurat <- subset(seurat, subset = DonorCategory != "Primary Infertility")

#setting celtype as a metadata column inorder not to loose it while setting donor category as ident
seurat@meta.data$celltypes <- Idents(seurat)
table(seurat@meta.data$celltypes)

#subset celltype for each celltype
celltype_seurat <- subset(seurat, subset = celltypes %in% c("Leydig "))
table(celltype_seurat@meta.data$DonorCategory)
Idents(celltype_seurat) <- celltype_seurat@meta.data$DonorCategory

# Extract metadata
meta <- celltype_seurat@meta.data

table(celltype_seurat@meta.data$Donor_Replicate)

# Get raw counts for DEseq analysis
counts <- GetAssayData(celltype_seurat, assay = "RNA",slot = "counts")
counts_t <- t(as.matrix(counts))

# Aggregate counts by sample_id
pseudobulk <- rowsum(as.matrix(counts_t), group = meta$Donor_Replicate)
pseudobulk_counts <- t(pseudobulk)


#coldata
patient <- c("Donor1_1","Donor1_2","Donor2_1","Donor2_2","Donor3_1","Donor3_2","I1_1","I1_2","I2_1","K1_1","K1_2")
condition <- c("Normal","Normal","Normal","Normal","Normal","Normal","Primary Infertility","Primary Infertility","Secondary Infertility","KS","KS")

#for use in the case of studying germ cells in secondary infertility
# patient <- c("Donor1_1","Donor1_2","Donor2_1","Donor2_2","Donor3_1","Donor3_2","I2_1")
# condition <- c("Normal","Normal","Normal","Normal","Normal","Normal","Secondary Infertility")


# Create coldata
coldata <- data.frame(
  row.names = colnames(pseudobulk_counts),
  condition = condition,
  patient = patient)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Normal")
dds <- DESeq(dds)

resultsNames(dds)

celltype <- "Leydig"
Condition_file <- "secinf_vs_Normal"
Condition_graph <- "Secondary Infertility Vs Normal"

# Order by adjusted p-value (FDR)
res <- results(dds, name = "condition_Secondary.Infertility_vs_Normal")
res <- res[order(res$padj), ]
res <- res[!is.na(res$padj), ]

res

# View the top differentially expressed genes
head(res)

rpgenes <- c("RPL39","RPL10","RPS2","RPL21","RPL18A","RPL13A","RPS17","RPS19","RPS18","RPS14","RPL34","RPS12","RPS27A","RPLP1","RPS8","RPL12","RPL23A","RPL18","RPS28","RPSA","RPS15A","RPS3","RPS21","RPL10A","RPL5","RPL24","RPL6","RPL9","RPL22","RPL14","RPL37","RPS16","RPL35A","RPL11","RPS7","RPS23","RPL26","RPL41","RPL7","RPL23","RPS20","RPL27","RPL31","RPL27A","RPS27","RPL3","RPLP0","RPS3A","RPS5","RPL4","RPS29","RPS11","RPL38","RPL37A","RPS9","RPL7A","RPL8","RPS6","RPL32","RPL13","RPS13","RPL36","RPLP2","RPS24","RPL15","RPL19","RPS15","RPL35","RPL30","RPS25","RPL28","RPS10","RPL29","RPS26","RPL17","RPL36A")


# Filter for significant genes (e.g., FDR < 0.05)
sig_res <- res[which(res$padj < 0.05), ]
res_filtered <- sig_res[!is.na(sig_res$padj) & sig_res$log2FoldChange>1 | sig_res$log2FoldChange<1, ]
filename <-paste0(celltype,"_",Condition_file,"_Deseq2_sigres.csv")
write.csv(res_filtered, file = filename)

# Filter res_filtered for genes in your rpgenes list
rp_filtered <- res_filtered[rownames(res_filtered) %in% rpgenes, ]

# Save to CSV
rp_filename <- paste0(celltype, "_", Condition_file, "_Deseq2_RPgenes.csv")
write.csv(rp_filtered, file = rp_filename)

res$gene <- rownames(res)
res$highlight <- ifelse(res$gene %in% rpgenes, "RP genes", "Other")

library(ggplot2)
library(ggrepel)

# Corrected regulation assignment
res$regulation <- ifelse(res$log2FoldChange > 1, "Upregulated",
                         ifelse(res$log2FoldChange < -1, "Downregulated", "Notsignificant"))

# New column for combined colors
res$combined_color <- ifelse(res$highlight == "RP genes",
                             ifelse(res$regulation == "Upregulated", "RP_Upregulated",
                                    ifelse(res$regulation == "Downregulated", "RP_Downregulated", "Insignificant")),
                             ifelse(res$regulation == "Upregulated", "Upregulated",
                                    ifelse(res$regulation == "Downregulated", "Downregulated", "Insignificant")))

# Color values
color_values <- c("RP_Upregulated" = "darkblue", "RP_Downregulated" = "red",
                  "Upregulated" = "lightskyblue3", "Downregulated" = "pink3", "Insignificant" = "grey")

# Shape value
shape_value <- c("RP genes" = 17, "Other" = 16)

plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = combined_color, shape = highlight), alpha = 0.7, size = 2) +  # Smaller size
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = shape_value) +
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed") +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(-15, 10)) +           # Set x-axis limits
  scale_y_continuous(limits = c(0, 90)) +             # Set y-axis limits
  theme_minimal(base_size = 14) +
  labs(
    title = paste0(celltype," ",Condition_graph," Volcano plot"),
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~p-value),
    color = "Regulation",
    shape = "Gene Group"
  ) +
  # Show gene names only for RP genes
  geom_text_repel(
    data = subset(res, highlight == "RP genes"),
    aes(label = gene),
    size = 3 , # Smaller size
    fontface = "bold", # Ensure it's not suppressed
    force = 4 
    ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
plot


file_name <- paste(celltype, "_Volcanoplot_",Condition_file,".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(plot)
dev.off()

file_path <- paste(celltype, "Volcanoplot",Condition_file,".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()


#PCA_to_see how the replicates are distributed

vsd <- vst(dds, blind = FALSE) 
pcaData <- plotPCA(vsd, intgroup = c("condition", "patient"), returnData = TRUE)
pcaData$patient <- colData(vsd)$patient
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = patient)) +
  geom_point(size = 3) +
  geom_text(aes(label = patient), nudge_x = 0.5, nudge_y = 0.5, size = 3) + # Add patient names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(paste0(celltype,"_PCA"))
pcaplot

file_path <- paste(celltype, "_normalVsdisease_PCA_plot",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(pcaplot)
# Close the PNG device
dev.off()

