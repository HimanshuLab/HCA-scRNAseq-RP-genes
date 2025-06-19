library(viridis)
library(ComplexHeatmap)
library(Seurat)
library(SeuratObject)
library(sp)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(EnhancedVolcano)
library(circlize)

rpgenes <- c("RPL23","RPS17","RPL36A","RPS20","RPL27A","RPL31","RPL7","RPL27","RPL4","RPS11","RPL38","RPS29","RPL13A","RPSA","RPL10A","RPL23A","RPS2","RPL21","RPS21","RPL18A","RPL18","RPL12","RPS8","RPS28","RPL41","RPL24","RPL9","RPS16","RPL22","RPL17","RPL14","RPS18","RPL6","RPS3","RPS9","RPL35","RPL3","RPS6","RPL5","RPS5","RPLP0","RPS10","RPS19","RPS25","RPS3A","RPS7","RPL26","RPS14","RPS23","RPL37","RPL36","RPL35A","RPS13","RPS27A","RPS12","RPL11","RPL34","RPS26","RPS15A","RPLP1","RPL37A","RPL15","RPL29","RPL7A","RPLP2","RPL8","RPS24","RPS27","RPL30","RPL32","RPS15","RPL19","RPL13","RPL28","RPL39","RPL10")

seurat <- readRDS("Trachea_14_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_trachea <- as.matrix(mean)
rownames(mean_trachea) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_trachea <- as.matrix(variance)
rownames(variance_trachea) <- rpgenes
IOD <- variance/mean
IOD_trachea <- as.matrix(IOD)
rownames(IOD_trachea) <- rpgenes

seurat <- readRDS("Stomach_05_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_stomach <- as.matrix(mean)
rownames(mean_stomach) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_stomach <- as.matrix(variance)
rownames(variance_stomach) <- rpgenes
IOD <- variance/mean
IOD_stomach <- as.matrix(IOD)
rownames(IOD_stomach) <- rpgenes

seurat <- readRDS("Spleen_06_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_spleen <- as.matrix(mean)
rownames(mean_spleen) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_spleen <- as.matrix(variance)
rownames(variance_spleen) <- rpgenes
IOD <- variance/mean
IOD_spleen <- as.matrix(IOD)
rownames(IOD_spleen) <- rpgenes

seurat <- readRDS("SmallIntestine_11_06.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_SI <- as.matrix(mean)
rownames(mean_SI) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_SI <- as.matrix(variance)
rownames(variance_SI) <- rpgenes
IOD <- variance/mean
IOD_SI <- as.matrix(IOD)
rownames(IOD_SI) <- rpgenes

seurat <- readRDS("Rectum_06_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_rectum <- as.matrix(mean)
rownames(mean_rectum) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_rectum <- as.matrix(variance)
rownames(variance_rectum) <- rpgenes
IOD <- variance/mean
IOD_rectum <- as.matrix(IOD)
rownames(IOD_rectum) <- rpgenes

seurat <- readRDS("Bile_15_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_bile <- as.matrix(mean)
rownames(mean_bile) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_bile <- as.matrix(variance)
rownames(variance_bile) <- rpgenes
IOD <- variance/mean
IOD_bile <- as.matrix(IOD)
rownames(IOD_bile) <- rpgenes

seurat <- readRDS("Bladder_16_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_bladder <- as.matrix(mean)
rownames(mean_bladder) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_bladder <- as.matrix(variance)
rownames(variance_bladder) <- rpgenes
IOD <- variance/mean
IOD_bladder <- as.matrix(IOD)
rownames(IOD_bladder) <- rpgenes

seurat <- readRDS("Muscle_05_04_recent.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_muscle <- as.matrix(mean)
rownames(mean_muscle) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_muscle <- as.matrix(variance)
rownames(variance_muscle) <- rpgenes
IOD <- variance/mean
IOD_muscle <- as.matrix(IOD)
rownames(IOD_muscle) <- rpgenes

seurat <- readRDS("Marrow_17_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_marrow <- as.matrix(mean)
rownames(mean_marrow) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_marrow <- as.matrix(variance)
rownames(variance_marrow) <- rpgenes
IOD <- variance/mean
IOD_marrow <- as.matrix(IOD)
rownames(IOD_marrow) <- rpgenes

seurat <- readRDS("Lymph_06_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_lymph <- as.matrix(mean)
rownames(mean_lymph) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_lymph <- as.matrix(variance)
rownames(variance_lymph) <- rpgenes
IOD <- variance/mean
IOD_lymph <- as.matrix(IOD)
rownames(IOD_lymph) <- rpgenes

seurat <- readRDS("esophagous_15_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_oesophagus <- as.matrix(mean)
rownames(mean_oesophagus) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_oesophagus <- as.matrix(variance)
rownames(variance_oesophagus) <- rpgenes
IOD <- variance/mean
IOD_oesophagus <- as.matrix(IOD)
rownames(IOD_oesophagus) <- rpgenes

seurat <- readRDS("Blood_19_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_blood <- as.matrix(mean)
rownames(mean_blood) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_blood <- as.matrix(variance)
rownames(variance_blood) <- rpgenes
IOD <- variance/mean
IOD_blood <- as.matrix(IOD)
rownames(IOD_blood) <- rpgenes

seurat <- readRDS("Skin_15_04.rds")
levels(seurat)
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_skin <- as.matrix(mean)
rownames(mean_skin) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_skin <- as.matrix(variance)
rownames(variance_skin) <- rpgenes
IOD <- variance/mean
IOD_skin <- as.matrix(IOD)
rownames(IOD_skin) <- rpgenes

seurat <- readRDS("Liver_12_04.rds")
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_liver <- as.matrix(mean)
rownames(mean_liver) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_liver <- as.matrix(variance)
rownames(variance_liver) <- rpgenes
IOD <- variance/mean
IOD_liver <- as.matrix(IOD)
rownames(IOD_liver) <- rpgenes


#for testis choosing one patient 
seurat <- readRDS("Testis_Adult_05_04.rds")
table(seurat@meta.data$Patient_Replicate)
seurat <- subset(seurat, subset = Patient_Replicate %in% c("Donor1_1","Donor1_2") )
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_testis1 <- as.matrix(mean)
rownames(mean_testis1) <- rpgenes
variance <- apply(mat, 1, var) %>% as.numeric()
variance_testis1 <- as.matrix(variance)
rownames(variance_testis1) <- rpgenes
IOD <- variance/mean
IOD_testis1 <- as.matrix(IOD)
rownames(IOD_testis1) <- rpgenes

colnames(IOD_blood) <- c("Blood")
colnames(IOD_oesophagus) <- c("Oesophagus")
colnames(IOD_liver) <- c("Liver")
colnames(IOD_lymph) <- c("Lymph")
colnames(IOD_marrow) <- c("Marrow")
colnames(IOD_muscle) <- c("Muscle")
colnames(IOD_bladder) <- c("Bladder")
colnames(IOD_bile) <- c("Bile")
colnames(IOD_rectum) <- c("Rectum")
colnames(IOD_SI) <- c("SmallIntestine")
colnames(IOD_skin) <- c("Skin")
colnames(IOD_stomach) <- c("Stomach")
colnames(IOD_testis1) <- c("Testis")
colnames(IOD_trachea) <- c("Trachea")
colnames(IOD_spleen) <- c("Spleen")


colnames(mean_blood) <- c("Blood")
colnames(mean_oesophagus) <- c("Oesophagus")
colnames(mean_liver) <- c("Liver")
colnames(mean_lymph) <- c("Lymph")
colnames(mean_marrow) <- c("Marrow")
colnames(mean_muscle) <- c("Muscle")
colnames(mean_bladder) <- c("Bladder")
colnames(mean_bile) <- c("Bile")
colnames(mean_rectum) <- c("Rectum")
colnames(mean_SI) <- c("SmallIntestine")
colnames(mean_skin) <- c("Skin")
colnames(mean_stomach) <- c("Stomach")
colnames(mean_testis1) <- c("Testis")
colnames(mean_trachea) <- c("Trachea")
colnames(mean_spleen) <- c("Spleen")


df_list_mean_data <- list(mean_blood,mean_oesophagus,mean_liver,mean_lymph,mean_marrow,mean_muscle,mean_bladder,mean_bile,mean_rectum,mean_SI,mean_skin,mean_stomach,mean_testis1,mean_trachea,mean_spleen) #mean_Skin_130973,
df_list_mean_data <- as.data.frame(df_list_mean_data)
df_list_mean <- as.matrix(df_list_mean_data)
max(df_list_mean)


# Define breaks for your data range (adjust as needed)
breaks <- seq(0, 4, length.out = 100)

# Create a color function using viridis palette
col_plasma <- colorRamp2(breaks, rev(viridis(length(breaks), option = "plasma")))

mean_ht <- Heatmap(df_list_mean_data, name = "mean",
                  cluster_columns = TRUE,
                  show_column_dend = TRUE,
                  cluster_row_slices = TRUE,
                  cluster_column_slices = TRUE,
                  column_title_gp = gpar(fontsize = 6),
                  column_gap = unit(0.5, "mm"),
                  cluster_rows = TRUE,
                  show_row_dend = TRUE,
                  col = col_plasma,
                  row_names_gp = gpar(fontsize = 7),
                  column_title_rot = 90,
                  show_column_names = TRUE,
                  use_raster = TRUE,
                  raster_quality = 4)

mean_ht


file_name <- paste("alltissue__mean_heatmap",".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(mean_ht)
dev.off()

file_path <- paste("alltissue_mean_heatmap",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(mean_ht)
# Close the PNG device
dev.off()

write.csv(df_list_mean_data,"mean_alltissues.csv")


df_list_IOD_data <- list(IOD_blood,IOD_oesophagus,IOD_liver,IOD_lymph,IOD_marrow,IOD_muscle,IOD_bladder,IOD_bile,IOD_rectum,IOD_SI,IOD_skin,IOD_stomach,IOD_testis1,IOD_trachea,IOD_spleen)#,IOD_skin2) #IOD_Skin_130973,
df_list_IOD_data <- as.data.frame(df_list_IOD_data)
df_list_IOD <- as.matrix(df_list_IOD_data)
write.csv(df_list_IOD_data,"IOD_alltissues.csv")


IOD_ht<- Heatmap(df_list_IOD_data, name = "IOD",
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        cluster_row_slices = TRUE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 6),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        col = colorRamp2(c(0,0.25,0.5,0.75,1,1.5,2), c("darkblue","lightskyblue","mistyrose","wheat","yellow","orange","darkred")),
        row_names_gp = gpar(fontsize = 7),
        column_title_rot = 90,
        #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))), 
        show_column_names = TRUE,
        use_raster = TRUE,
        raster_quality = 4)

IOD_ht
file_name <- paste("alltissue_heatmap_IOD",".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(IOD_ht)
dev.off()

file_path <- paste("alltissue_heatmap_IOD",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(IOD_ht)
# Close the PNG device
dev.off()


#Defining the percentiles and finding the high IOD values
hist_data <- unlist(df_list_IOD_data)

# Calculate percentiles
p10 <- quantile(hist_data, 0.10)
p25 <- quantile(hist_data, 0.25)
p50 <- quantile(hist_data, 0.50)
p75 <- quantile(hist_data, 0.75)
p90 <- quantile(hist_data, 0.90)

png("IOD_density_plot.png", width = 1200, height = 800, res = 150)
# Plot the density with a black main line and custom axis labels
plot(density(hist_data), 
     main = "Distribution of Gene Expression IOD Values", 
     xlab = "Gene Expression (IOD)", 
     ylab = "Density", 
     col = "black", 
     lwd = 2)

# Add percentile lines
abline(v = p10, col = "darkblue", lwd = 2, lty = 1)
abline(v = p25, col = "grey40", lwd = 2, lty = 1)
abline(v = p50, col = "red", lwd = 2, lty = 1)
abline(v = p75, col = "grey40", lwd = 2, lty = 1)
abline(v = p90, col = "darkblue", lwd = 2, lty = 1)

legend("topright", 
       legend = c("10th percentile", "25th percentile", "Median", "75th percentile", "90th percentile"),
       col = c("darkblue", "grey40", "red", "grey40", "darkblue"),
       lty = c(1, 1, 1, 1, 1),
       lwd = 2,
       bty = "n",
       y.intersp = 0.7)
dev.off()

svg("IOD_density_plot.svg", width = 10, height = 6)
# Plot the density with a black main line and custom axis labels
plot(density(hist_data), 
     main = "Distribution of Gene Expression IOD Values", 
     xlab = "Gene Expression (IOD)", 
     ylab = "Density", 
     col = "black", 
     lwd = 2)

# Add percentile lines
abline(v = p10, col = "darkblue", lwd = 2, lty = 1)
abline(v = p25, col = "grey40", lwd = 2, lty = 1)
abline(v = p50, col = "red", lwd = 2, lty = 1)
abline(v = p75, col = "grey40", lwd = 2, lty = 1)
abline(v = p90, col = "darkblue", lwd = 2, lty = 1)

legend("topright", 
       legend = c("10th percentile", "25th percentile", "Median", "75th percentile", "90th percentile"),
       col = c("darkblue", "grey40", "red", "grey40", "darkblue"),
       lty = c(1, 1, 1, 1, 1),
       lwd = 2,
       bty = "n",
       y.intersp = 0.2)
dev.off()

 
# Correctly count number of entries above column-wise 90th percentile
counts_above_p90 <- sapply(1:ncol(df_list_IOD_data), function(i) {
  sum(df_list_IOD_data[, i] > p90)
})

# Name the result properly
names(counts_above_p90) <- colnames(df_list_IOD_data)
counts_sorted <- sort(counts_above_p90, decreasing = TRUE)

# Choose a colorblind-friendly palette (1 color per bar)
colors <- viridis(length(counts_sorted), option = "D")

# Create the barplot
png("RP_density_tissuewise_IOD.png", width = 1200, height = 800, res = 150)
bp <- barplot(counts_sorted,
              main = "Entries above 90th percentile per tissue",
              ylab = "Count",
              col = colors,
              xaxt = "n",
              ylim = c(0, max(counts_sorted) + 5))  # add space for labels

# Add rotated x-axis labels
text(x = bp,
     y = par("usr")[3] - 0.5,
     labels = names(counts_sorted),
     srt = 45,
     adj = 1,
     xpd = TRUE,
     cex = 0.8)

# Add count labels above each bar
text(x = bp,
     y = counts_sorted + 1,
     labels = counts_sorted,
     cex = 0.8,
     font = 2)
dev.off()

svg("RP_density_tissuewise_IOD.svg", width = 10, height = 6)
bp <- barplot(counts_sorted,
              main = "Entries above 90th percentile per tissue",
              ylab = "Count",
              col = colors,
              xaxt = "n",
              ylim = c(0, max(counts_sorted) + 5))  # add space for labels

# Add rotated x-axis labels
text(x = bp,
     y = par("usr")[3] - 0.5,
     labels = names(counts_sorted),
     srt = 45,
     adj = 1,
     xpd = TRUE,
     cex = 0.8)

# Add count labels above each bar
text(x = bp,
     y = counts_sorted + 1,
     labels = counts_sorted,
     cex = 0.8,
     font = 2)
dev.off()

#heatmap of different patient replicates of testis

seurat <- readRDS("Testis_Adult_05_04.rds")
table(seurat@meta.data$Patient_Replicate)
P1 <- subset(seurat, subset = Patient_Replicate %in% c("Donor1_1","Donor1_2") )
mat <- P1[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_testis1 <- as.matrix(mean)
rownames(mean_testis1) <- rpgenes
colnames(mean_testis1) <- "Donor_1"
variance <- apply(mat, 1, var) %>% as.numeric()
variance_testis1 <- as.matrix(variance)
rownames(variance_testis1) <- rpgenes
IOD <- variance/mean
IOD_testis1 <- as.matrix(IOD)
rownames(IOD_testis1) <- rpgenes
colnames(IOD_testis1) <- "Donor_1"

P2 <- subset(seurat, subset = Patient_Replicate %in% c("Donor2_1","Donor2_2") )
mat <- P2[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_testis2 <- as.matrix(mean)
rownames(mean_testis2) <- rpgenes
colnames(mean_testis2) <- "Donor_2"
variance <- apply(mat, 1, var) %>% as.numeric()
variance_testis2 <- as.matrix(variance)
rownames(variance_testis2) <- rpgenes
IOD <- variance/mean
IOD_testis2 <- as.matrix(IOD)
rownames(IOD_testis2) <- rpgenes
colnames(IOD_testis2) <- "Donor_2"

P3 <- subset(seurat, subset = Patient_Replicate %in% c("Donor3_1","Donor3_2") )
mat <- P3[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
mean_testis3 <- as.matrix(mean)
rownames(mean_testis3) <- rpgenes
colnames(mean_testis3) <- "Donor_3"
variance <- apply(mat, 1, var) %>% as.numeric()
variance_testis3 <- as.matrix(variance)
rownames(variance_testis3) <- rpgenes
IOD <- variance/mean
IOD_testis3 <- as.matrix(IOD)
rownames(IOD_testis3) <- rpgenes
colnames(IOD_testis3) <- "Donor_3"



df_list_mean_testis <- list(mean_testis1,mean_testis2,mean_testis3) #mean_Skin_130973,
df_list_mean_testis <- as.data.frame(df_list_mean_testis)
df_list_mean_testis <- as.matrix(df_list_mean_testis)

write.csv(df_list_mean_testis,"mean_testis.csv")



# Define breaks for your data range (adjust as needed)
breaks <- seq(0, 4, length.out = 100)

# Create a color function using viridis palette
col_plasma <- colorRamp2(breaks, rev(viridis(length(breaks), option = "plasma")))

mean_ht <- Heatmap(df_list_mean_testis, name = "mean",
                   cluster_columns = FALSE,
                   show_column_dend = TRUE,
                   cluster_row_slices = TRUE,
                   cluster_column_slices = TRUE,
                   column_title_gp = gpar(fontsize = 6),
                   column_gap = unit(0.5, "mm"),
                   cluster_rows = TRUE,
                   show_row_dend = TRUE,
                   col = col_plasma,
                   row_names_gp = gpar(fontsize = 7),
                   column_title_rot = 90,
                   show_column_names = TRUE,
                   use_raster = TRUE,
                   raster_quality = 4,
                   heatmap_width = unit(40, "mm"))

mean_ht

df_list_IOD_testis <- list(IOD_testis1,IOD_testis2,IOD_testis3)
df_list_IOD_testis <- as.data.frame(df_list_IOD_testis)
df_list_IOD_testis <- as.matrix(df_list_IOD_testis)
write.csv(df_list_IOD_testis,"IOD_testis.csv")

IOD_ht<- Heatmap(df_list_IOD_testis, name = "IOD",
                 cluster_columns = FALSE,
                 show_column_dend = TRUE,
                 cluster_row_slices = TRUE,
                 cluster_column_slices = TRUE,
                 column_title_gp = gpar(fontsize = 6),
                 column_gap = unit(0.5, "mm"),
                 cluster_rows = FALSE,
                 show_row_dend = TRUE,
                 col = colorRamp2(c(0,0.25,0.5,0.75,1,1.5,2), c("darkblue","lightskyblue","mistyrose","wheat","yellow","orange","darkred")),
                 row_names_gp = gpar(fontsize = 7),
                 column_title_rot = 90,
                 #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))), 
                 show_column_names = TRUE,
                 use_raster = TRUE,
                 raster_quality = 4,
                 heatmap_width = unit(25, "mm"))

IOD_ht
file_name <- paste("testis_heatmap_IOD",".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(IOD_ht)
dev.off()

file_path <- paste("Testis_heatmap_IOD",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(IOD_ht)
# Close the PNG device
dev.off()# Correctly count number of entries above column-wise 90th percentile


counts_above_p90 <- sapply(1:ncol(df_list_IOD_testis), function(i) {
  sum(df_list_IOD_testis[, i] > p90)
})

# Name the result properly
names(counts_above_p90) <- colnames(df_list_IOD_testis)

# Choose a colorblind-friendly palette (1 color per bar)
colors <- viridis(length(counts_above_p90), option = "D")

# Create the barplot
png("RP_density_IOD_testis.png", width = 1200, height = 800, res = 150)
bp <- barplot(counts_above_p90,
              main = "Entries above 90th percentile per tissue",
              ylab = "Count",
              col = colors,
              xaxt = "n",
              width = 0.1,    # relative bar width (try 0.3-0.5)
              space = 0.1,    # reduce space between bars
              ylim = c(0, max(counts_above_p90) + 5))

# Add rotated x-axis labels
text(x = bp,
     y = par("usr")[3] - 0.5,
     labels = names(counts_above_p90),
     srt = 45,
     adj = 1,
     xpd = TRUE,
     cex = 0.8)

# Add count labels above each bar
text(x = bp,
     y = counts_above_p90 + 1,
     labels = counts_above_p90,
     cex = 0.8,
     font = 2)
dev.off()

svg("RP_density_testis_IOD.svg", width = 10, height = 6)
bp <- barplot(counts_above_p90,
              main = "Entries above 90th percentile per tissue",
              ylab = "Count",
              col = colors,
              width = 0.5,
              xaxt = "n",
              ylim = c(0, max(counts_above_p90) + 5))  # add space for labels

# Add rotated x-axis labels
text(x = bp,
     y = par("usr")[3] - 0.5,
     labels = names(counts_above_p90),
     srt = 45,
     adj = 1,
     xpd = TRUE,
     cex = 0.8)

# Add count labels above each bar
text(x = bp,
     y = counts_above_p90 + 1,
     labels = counts_above_p90,
     cex = 0.8,
     font = 2)
dev.off()

C1 <- as.data.frame(table(Idents(P1)))
colnames(C1) <- c("Cell types","P1")
C2 <- as.data.frame(table(Idents(P2)))
colnames(C2) <- c("Cell types","P2")
C3 <- as.data.frame(table(Idents(P3)))
colnames(C3) <- c("Cell types","P3")

testis_freq <- cbind(C1,C2[,2],C3[,2])
colnames(testis_freq) <- c("Cell types","P1","P2","P3")

library(ggplot2)
library(tidyr)
library(dplyr)


cell_data_long <- testis_freq %>%
  pivot_longer(cols = starts_with("P"), names_to = "Patient", values_to = "Count") %>%
  group_by(`Cell types`) %>%
  mutate(Total = sum(Count),
         Percentage = Count / Total,
         PercentLabel = paste0(round(Percentage * 100), "%"),
         Cell_label = paste0(`Cell types`, "; n=", Total)) %>%
  ungroup()

# Custom colors for each patient


# Plot with percentage labels inside bars
plot <- ggplot(cell_data_long, aes(x = reorder(Cell_label, Total), y = Percentage, fill = Patient)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  geom_text(aes(label = PercentLabel),
            position = position_fill(vjust = 0.5),
            color = "white", size = 5, fontface = "bold") +
  scale_fill_manual(values = my_colors,
                    labels = c("P1" = "Donor_1", "P2" = "Donor_2", "P3" = "Donor_3")) +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1.0)) +
  labs(x = "", y = "Percentage", fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80"),   # Optional: make gridlines black
    panel.grid.minor = element_line(color = "grey80"),   # Optional: make gridlines black
    axis.line = element_line(color = "black"),          # Optional: make axis lines black
    panel.border = element_rect(color = "black", fill = NA) # Optional: add a black border
  ) +
  coord_flip()

plot

file_name <- paste("testis_donor_celltypes",".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(plot)
dev.off()

file_path <- paste("testis_donor_celltypes",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(plot)
# Close the PNG device
dev.off()
plot
