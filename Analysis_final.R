library(ComplexHeatmap)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(EnhancedVolcano)
library(circlize)
library(viridis)
library(car)


filename <-"SmallIntestine_11_06.rds"
Tissue_name <- "SI"

seurat <- readRDS(filename)
rpgenes <- c("RPL23","RPS17","RPL36A","RPS20","RPL27A","RPL31","RPL7","RPL27","RPL4","RPS11","RPL38","RPS29","RPL13A","RPSA","RPL10A","RPL23A","RPS2","RPL21","RPS21","RPL18A","RPL18","RPL12","RPS8","RPS28","RPL41","RPL24","RPL9","RPS16","RPL22","RPL17","RPL14","RPS18","RPL6","RPS3","RPS9","RPL35","RPL3","RPS6","RPL5","RPS5","RPLP0","RPS10","RPS19","RPS25","RPS3A","RPS7","RPL26","RPS14","RPS23","RPL37","RPL36","RPL35A","RPS13","RPS27A","RPS12","RPL11","RPL34","RPS26","RPS15A","RPLP1","RPL37A","RPL15","RPL29","RPL7A","RPLP2","RPL8","RPS24","RPS27","RPL30","RPL32","RPS15","RPL19","RPL13","RPL28","RPL39","RPL10")

#Analysis 1
#Heterogenity of RP genes within a tissue using heatmap
#Extracting the expression matrix form seurat object
exp_matrix<- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
expvar <- apply(exp_matrix,1,var)
mean <- apply(exp_matrix,1,mean)
IOD <- expvar/mean
min(exp_matrix)
max(exp_matrix)
min(IOD)
max(IOD)

library(circlize)
# These limits of colours are common across tissues
col_list_expression = colorRamp2(c(0,1,2,3,4,5,6),c("#ffffff","#fff47f","#edb642","#8b6000","#522d00","#3b1500","#000000"))
col_list_variance <- list(variance = colorRamp2(c(0, 0.4, 0.8, 1.2,1.6,2.4,2.8), 
                                                c("darkblue", "lightskyblue","pink", "mistyrose", "yellow","orange", "darkred")))
col_list_IOD <- list(IOD = colorRamp2(c(0,0.25,0.5,0.75,1,1.5,2), c("darkblue","lightskyblue","mistyrose","wheat","yellow","orange","darkred")))


ha <- HeatmapAnnotation(IOD = IOD, which=c("row"),title = NULL, col = col_list_IOD)
ht <- Heatmap(exp_matrix, name = paste("RP Expression"),  
              column_split = Idents(seurat), 
              cluster_columns = TRUE,
              show_column_dend = TRUE,
              cluster_row_slices = TRUE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontSize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = FALSE,
              show_row_dend = TRUE,
              column_title_rot = 90,
              col = col_list_expression,
              row_names_gp = gpar(fontsize = 5),
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))), 
              left_annotation = ha,
              show_column_names = FALSE,
              use_raster = TRUE,
              raster_quality = 4)
ht

file_name <- paste(Tissue_name,"heatmap",".svg", sep = "")
svg(file_name,width = 10, height = 8)
print(ht)
dev.off()

file_path <- paste(Tissue_name,"heatmap",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
print(ht)
# Close the PNG device
dev.off()

#Analysis 2
#Differentially expressed genes across patients in the tissues studied when there are multiple patients in the dataset
levels(seurat)

#for tissues with multiple donors use this
#----
seurat <- PrepSCTFindMarkers(seurat)


# Defining the celltypes that are to be compared
#celltype_list <- c("SSC","DS","Early PS","Late PS","Round Sptd","Elongated Sptd","Sperm") 
celltype_list <- c("Endothelial","Leydig","Myoid","Macrophage")


# Create output directory if it doesn't exist
if (!dir.exists("conserved_markers_Infant_testis")) dir.create("conserved_markers_Infant_testis")
if (!dir.exists("conserved_markers_Infant_testis/RP_markers")) dir.create("conserved_markers_Infant_testis/RP_markers")

#DEG analysis for dataest with multiple donors
celltype_list <- levels(seurat)
levels(seurat)
seurat@meta.data$celltypes <- Idents(seurat)

for (celltype in celltype_list) {
  # Skip if celltype is the only one left
  if (length(celltype_list) < 2) next 
  # Defining ident_2 as all cell types except the celltype compared
  ident_2 <- setdiff(celltype_list, celltype)
  # Find conserved markers
  markers <- FindConservedMarkers(
    seurat,
    ident.1 = celltype,
    ident.2 = ident_2,
    grouping.var = "Patient_Replicate", #the metadata column where the donor information is stored
    assay = "SCT",
    slot = "data",
    min.cells.group = 3,
    meta.method = metap::minimump,
    verbose = TRUE
  )
  # extracting the p values from markers and corrceting them for multipl testing
  min_pvals <- markers$minimump_p_val
  
  # Apply Bonferroni correction
  adjusted_pvals <- p.adjust(min_pvals, method = "bonferroni")
  
  # Combine adjusted p-values with the original data
  markers$adjusted_pval <- adjusted_pvals
  
  markers$minimum_avg_foldchange <- apply(markers[, grep("avg_log2FC", colnames(markers))], 1, min)
  markers$maximum_avg_foldchange <- apply(markers[, grep("avg_log2FC", colnames(markers))], 1, max)
  
  markers$gene <- rownames(markers)
  graphname <- paste0("Volcano plots of ",celltype ," VS other celltypes")
  filename <- paste0("conserved_markers_Infant_testis/",Tissue_name,"_",celltype, "_vs_others_unfiltered.csv")
  write.csv(markers, filename)
  message(paste("Saved subsetted markers to:", filename))
  
  
  # Skip if no markers found
  if (nrow(markers) == 0) next
  
  #not applying the condition Pct1 > Pct2 as only the upregulated gene will show
  
  plot <- EnhancedVolcano(markers, 
                          rownames(markers),
                          x ="minimum_avg_foldchange",
                          y ="adjusted_pval",
                          selectLab = rpgenes,
                          title = graphname,
                          subtitle = 'Highlighting RP genes')
  
  
  plot
  
  # choosing those genes for which the signs of all foldchanges are similar
  significant_filtered_markers <- markers %>%
    filter(
      adjusted_pval < 0.05,
      minimum_avg_foldchange > 1 & maximum_avg_foldchange > 1 |
        minimum_avg_foldchange < -1 & maximum_avg_foldchange < -1 
    ) 
  
  # Filter the markers based on the gene list
  subset_markers <- significant_filtered_markers[rownames(significant_filtered_markers) %in% rpgenes, ]
  
  # Check if any genes exist in the markers
  if (nrow(subset_markers) == 0) {
    print(celltype)
    message("No RP genes from the list exist in the markers.")
  } else {
    # Save the subsetted markers to a CSV file
    filename <- paste0("conserved_markers_Infant_testis/RP_markers/",Tissue_name,"_",celltype, "_vs_others_RP_rpsubset.csv")
    write.csv(subset_markers, filename)
    message(paste("Saved subsetted markers to:", filename))
  }
  
  # Save results if markers remain after filtering
  if (nrow(significant_filtered_markers) > 0) {
    filename <- paste0("conserved_markers_Infant_testis/", Tissue_name,"_",celltype, "_vs_others_significant_markers.csv")
    write.csv(significant_filtered_markers, filename)
    message("Saved: ", filename)
  }
}
#----

# finding DEGs for tissues with single patient

#----
markers <- FindAllMarkers(
  seurat,
  assay = "SCT",
  slot = "data",
  verbose = TRUE
)

significant_filtered_markers <- markers %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 1,  # Better practice than separate >1/< -1
    !is.na(p_val_adj)     # Exclude NA values if any
  )

# Filter the markers based on the gene list
subset_markers <- significant_filtered_markers[significant_filtered_markers$gene %in% rpgenes, ]

filename <- paste0(Tissue_name,"_significant_markers.csv")
write.csv(significant_filtered_markers, filename)

filename <- paste0(Tissue_name,"_RP_markers.csv")
write.csv(subset_markers, filename)
#----

#comparing RP genes with ubiquitously expressed genes

all.genes <- rownames(seurat)

dna_dsbreak_repair <- c("GEN1","HUS1","INIP","RECQL","LIG3","ZFYVE26","MORF4L1","APLF","MRE11","UBE2N","DCLRE1C","DCLRE1B","ATP23","PNKP","RNF8","YY1","SLX4","NSD2","NBN","SMARCAD1","NHEJ1","DCLRE1A","HELQ","XRCC5","XRCC6","TOP3A","INTS3","XRCC4","INO80","RNF169","RNF168","SEM1","SETX","TDP2","NABP2","WDR48","VCP","PRKDC","SMC6","SMC5","MBTD1","EPC1","TP53BP1","ACTR8","NPC1","SPIDR","POLA1","RAD51D","SFPQ","ATM","SWSAP1","PSMD14","SWI5","HMGB1","CHD4","HMGB2","MRPL36","BABAM2","RAD21","BABAM1","UIMC1","CYREN","KDM2A","SLF2","EYA3","POLN","SLF1","NSMCE3","PARP1","NSMCE2","RPA2","NSMCE1","RPA3","RPA1","RAD50","ERCC1","RAD52","ERCC4","AUNIP","UBE2V2","ZMYND8","MCM9","NUCKS1","POLB","WRN","ZSWIM7","REV3L") #Double-Strand Break Repair (GO:0006302)
common_dnarepair <- intersect(dna_dsbreak_repair,all.genes)

cellular_resp <- c("NDUFC2","NDUFC1","SDHB","COX6C","SDHC","SDHA","NDUFS8","LYRM7","NDUFS6","UQCRQ","NDUFS5","NDUFS4","UQCRC2","NDUFS3","UQCRC1","CYCS","NDUFS2","NDUFS1","COX7C","NDUFB9","COX7B","NDUFB8","NDUFB11","NDUFB7","NDUFB6","NDUFB5","NDUFB10","NDUFB4","NDUFB3","UQCR10","ETFDH","NDUFB2","NDUFB1","UQCR11","COX6A1","UQCRFS1","COX4I1","ETFRF1","NDUFV3","NDUFV2","NDUFV1","COX8A","NDUFA9","NDUFA8","NDUFA7","SURF1","MDH2","NDUFA6","NDUFA5","NDUFA4","NDUFA3","NDUFA2","NDUFA1","COX6B1","OXA1L","NDUFAB1","COX5A","NDUFA12","NDUFA13","UQCRB","NDUFA10","NDUFA11","MTFR1L","ETFA","ETFB","COX5B","COX7A2","UQCRH","CHCHD5","BLOC1S1","COX15") #Cellular Respiration (GO:0045333)
common_cellular <- intersect(cellular_resp,all.genes)

endosomal_transport <- c("ANKRD27","PLEKHA3","ANKFY1","ZFYVE16","ACAP2","UBE2O","TBC1D23","SQSTM1","VPS29","TMEM87A","SNX27","VTI1B","STX6","VTI1A","PLEKHJ1","REPS1","SNF8","PICALM","CHMP2B","CHMP3","CHMP2A","VAMP3","WASHC5","DPY30","ARHGAP1","ARL1","GCC2","VPS26C","VPS26A","CORO1C","BLOC1S1","BLOC1S2","CHMP1B","CCDC93","WASHC3","WASHC4","RAB6D","EIPR1","CMTM6","RBSN","HEATR5B","PIKFYVE","HEATR5A","RAB7A","VPS35L","VPS4A","DENND1B","VPS4B","CLTC","AP5M1","SNX4","SNX5","SNX2","VPS50","SNX3","DENND10","VPS52","SNX1","VPS53","EPS15","EVI5","RAB6C","SNX6","RAB6A","SPAG9","RAB10","SNX17","TBC1D5","RAB14","LAMTOR1","RAB9A","ARF6","GOSR1","ALMS1","DCTN1","CLN5","ERC1","TRIM27","VPS35") #Endosomal Transport (GO:0016197)
common_endosomal <- intersect(endosomal_transport,all.genes)

ERtogolgi <- c("SEC24C","RAB2A","SEC24B","YIPF4","TRAPPC11","YIPF5","TRAPPC10","NRBP1","YIPF7","MPPE1","CNIH4","RAB1A","PGAP1","KDELR1","KDELR2","STX17","PDCD6","TFG","TRAPPC13","VTI1A","TRAPPC2L","TRAPPC3","TRAPPC4","TRAPPC2","TRAPPC8","TRAPPC2B","HYOU1","PEF1","COPA","COPB2","COPB1","CUL3","KLHL12","SEC31A","SEC13","SEC16A","COPG1","ERGIC3","ERGIC2","ARF4","TMED10","SAR1A","SAR1B","WHAMM","IER3IP1","SPAST","SCFD1","TMED4","TMED3","TMED2","TMED7","TMED6","TMED5","SEC22A","BCAP31","ARCN1","TRAPPC6B","VAPA","VAPB","TRIP11","USO1","GOPC","SEC22C","SEC22B","BET1","GOSR1","ATL2","ATL3","GOSR2","MIA2","MIA3","USE1","LMAN1","LMAN2","BCAP29") #Endoplasmic Reticulum To Golgi Vesicle-Mediated Transport (GO:0006888)
common_ERtogolgi <- intersect(ERtogolgi,all.genes)

rnamet <- c("DICER1","HNRNPK","DDX39B","HNRNPL","DKC1","XRN2","HNRNPF","HNRNPD","ZNF830","RNASEH2C","RNASEH2B","HNRNPU","NOL9","PARN","TBRG4","EXOSC8","EXOSC7","EXOSC10","EXOSC6","DHX36","EXOSC9","PABPN1","EXOSC3","EXOSC2","EXOSC1","DNTTIP2","TFIP11","RNASEH1","SUPV3L1","SETX","PRPF4","FTO","DIS3L","MPHOSPH10","SYNCRIP","PCBP2","ALKBH1","PAPOLG","PAPOLA","PNPT1","DIS3","NCBP1","AGO2","AGO3","AGO1","PABPC4","DDX3X","NUFIP1","DHX8","FASTKD1","DHX9","FASTKD2","BICD1","ADAR","ATXN2","DUSP11","RBM39","RNGTT","DDX17","THOC1","HNRNPUL1","SNRNP40","HNRNPH1","HNRNPH3","DDX24","NUDT3","DDX20","RBM4","RBM5","RBM3","PPP1R8","MTREX","TENT2","RBM6") #RNA Metabolic Process (GO:0016070)
common_rnamet <- intersect(rnamet,all.genes)

rnapoltrans <- c("MED1","MED4","SUPT16H","DEK","GTF2F1","GTF2F2","MED7","XRN2","COPS2","SSU72","SUB1","TBPL1","CDK12","CCNK","CCNT2","CCNT1","CCNH","GTF2B","GTF2E2","SCAF8","TAF1A","E2F3","LEO1","RPRD1B","PKNOX1","XBP1","RPRD1A","ZNF141","PCID2","MED27","SUPT5H","MED23","CREB1","TFDP1","SETX","MED20","MNAT1","ATF4","ICE2","CBFB","ICE1","CTR9","NCBP2","NCBP3","NCBP1","PEX2","CDK7","NFX1","GTF2A1","GTF2A2","ELOF1","NFAT5","PMF1","SUPT4H1","RTF1","ZC3H8","MED17","CDC73","PCF11","SNAPC3","SNAPC5","SNAPC1","ZNF345","PAAF1","TBP","TAF13","PBRM1","ADRM1","TAF11","PARP1","TAF12","ELP2","GTF2H1","TAF10","ELP4","GTF2H2","GTF2H3","GTF2H5","ERCC3","TAF9","TRIP11","TAF8","ERCC6","TAF7","TCEA1","RBMX","TAF5","TAF3","TAF2","DDX21","POLR2B","POLR2C","POLR2D","POLR2E","ELOA","POLR2F","ELOB","POLR2G","ELOC","POLR2H","POLR2I","POLR2J","POLR2K","SETD2","POLR2L","TRIM24")#Transcription By RNA Polymerase II (GO:0006366)
common_rnapoltrans <- intersect(rnapoltrans,all.genes)

phosporylation <- c("PRKAB1","ATG14","ATG13","CDK16","TESK2","PRKAA1","PPP4R1","LATS1","GMFB","CTBP1","PRKDC","TBK1","XYLB","CSNK1A1","BRAF","GPI","CDKL3","BRD2","CSNK2A2","CSNK2A1","STAT3","ADAM10","TAOK1","PKN2","CDKL1","SMG1","PRKAG1","RICTOR","JAK1","PRPF4B","VRK2","HIPK1","ILF3","COPS2","COPS8","DGKH","DGKE","ROCK1","ROCK2","CCNT1","MAP2K1","DYRK1A","STK24","MAK","BIRC6","MAP3K12","MAP3K13","NEK1","NEK4","EIF2AK2","EIF2AK4","PIK3CB","PIK3CA","CDK4","PIK3C3","TOP1","BMPR2","LRRK1","MAPKAP1","FGGY","PGK1","CHUK","MTOR","BECN1","TLK2","PANK2","CSNK1E","SRPK1","SRPK2","WNK1","CSNK2B","PDK1","PDK3","ATM","FAM20B","ATR","STK38","MAP4K5","AURKC","PAK2","MORC3","CREB1","TBCK","DYRK3","DYRK2","MARK3","ACVR2A","CLK3","FER","GALK2","RSKR","GSK3B","ILK","PIK3R4","STK3","STK4","GRK4","MKNK1","TBP","CDK11B","CDK11A","TAF13","PDPK1","TAF11","TAF12","TAF10","OXSR1","YWHAZ","PEAK1","TAF9","TAF8","TAF7","TAF5","TAF3","TAF2","RB1CC1","RPS6KA5","GK5","MAPK9","MAPK1") #Phosphorylation (GO:0016310)
common_phosporylation <- intersect(phosporylation,all.genes)

chromatin_remodeling <- c("PIH1D1","DEK","GATAD2A","SMARCA4","GATAD2B","SMARCA5","APLF","ASF1A","HNRNPC","CGGBP1","RNF2","ZNF274","IWS1","RNF8","MTA3","SF3B1","PHC1","DAXX","ACTL6A","INO80","BAZ1B","HCFC2","SRPK2","UCHL5","PAXIP1","BAZ2A","ZNHIT1","DPF2","MIER1","SKP1","SS18","MCRS1","PHF10","TADA2A","MTF2","CTCF","CHRAC1","BABAM1","SUDS3","LIN54","NPM1","PBRM1","CHD1L","PER2","RYBP","TRIP12","WAC","HDGF","SETD4","RSF1","NUDT5","PCGF1","PWWP2A","NFRKB","BRD7","SETD2","USP3","POLE3","CHTOP","KDM1A","MACROH2A1","RUVBL1","YY1","RBBP4","RUVBL2","RBBP7","SMARCAD1","MORC2","BPTF","GLMN","PCID2","RNF168","INO80D","INO80C","GATAD1","ARID2","KDM4A","TTI1","DR1","ACTR8","SETDB2","SETDB1","CBX3","SMARCE1","SFPQ","TOP1","KDM5A","HP1BP3","NFAT5","CHD8","PSIP1","UBR5","CHD6","CHD4","CHD2","CHD1","SMARCC1","SMARCC2","ATRX","ARID1A","ARID1B","CFDP1","MCM3AP","REST","ERCC6","USP15","HDAC4","HDAC2","DDX23","DDX21","PRIMPOL","SART3","TDG","GET1") #Chromatin Remodeling (GO:0006338)
common_chromatin_remodeling <- intersect(chromatin_remodeling,all.genes)

genepathway_list <- list(common_dnarepair,common_cellular,common_chromatin_remodeling,common_endosomal,common_ERtogolgi,common_rnamet,common_rnapoltrans,common_phosporylation)

#Identifying IOD of RP genes throughout then tissue
mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
variance <- apply(mat, 1, var) %>% as.numeric()
final_variance_matrix_UEG <- as.matrix(variance)
rownames(final_variance_matrix_UEG) <- rpgenes
var_mat_melt <- melt(final_variance_matrix_UEG)
var_mat_melt$Var2 <- "RP genes"
IOD <- variance/mean
IOD_matrix <- as.matrix(IOD)
finalIOD_matrix_UEG <- IOD_matrix
rownames(finalIOD_matrix_UEG) <- rpgenes
IOD_mat_melt <- melt(finalIOD_matrix_UEG)
IOD_mat_melt$Var2 <- "RP genes"

#finding IOD of all the other genes and adding them to IOD_mat_melt
list <- c("Dna_repair","Cellular_respiration","Chromatin_remodeling","Endosomal_transport","ERtogolgi","RNA_met","RNA_poltrans","Phosporylation")
for (j in 1:length(genepathway_list)) {
  
  genelist <- as.list(genepathway_list[j])
  genelist <- unlist(genelist)
  mat <- seurat[["SCT"]]$data[genelist, ] %>% as.matrix()
  mean <- apply(mat, 1, mean) %>% as.numeric()
  variance <- apply(mat, 1, var) %>% as.numeric()
  var_matrix <- as.matrix(variance)
  rownames(var_matrix) <- genelist
  var_mat_melt_j <- melt(var_matrix)
  var_mat_melt_j$Var2 <- list[j]
  IOD <- variance/mean
  IOD_matrix <- as.matrix(IOD)
  rownames(IOD_matrix) <- genelist
  IOD_mat_melt_j <- melt(IOD_matrix)
  IOD_mat_melt_j$Var2 <- list[j]
  
  var_mat_melt <- rbind(var_mat_melt,var_mat_melt_j)
  IOD_mat_melt <- rbind(IOD_mat_melt,IOD_mat_melt_j)
}

colnames(IOD_mat_melt) <- c("Genes","Pathways","IOD")
colnames(var_mat_melt) <- c("Genes","Pathways","variance")



plot_UEG <- ggplot(IOD_mat_melt, 
                   aes(x = Pathways, y = IOD, fill = Pathways)) + 
  geom_violin() +
  ylab("Index of Dispersion") +
  xlab("Genes") +
  ggtitle("Comparison of Index of dispersion of RP genes VS Ubiquitous gene sets") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d(option = "A")  # Discrete color mapping

plot_UEG

file_path <- paste(Tissue_name, "RP Vs UEG",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(plot_UEG)
# Close the PNG device
dev.off()

#levene's test to test overall difference
levene_result <- leveneTest(IOD ~ Pathways, data = IOD_mat_melt)

# Capture the result as a character string
result_text <- paste(capture.output(levene_result), collapse = "\n")

# Write the result to a text file with a message
write(paste("\nLevene Test Result RP Vs UEG:\n", result_text), 
      file = paste0(Tissue_name,"_","levene_test_result.txt"), 
      append = FALSE)

#Filtering the data to remove RPs and levene's test to see comparisions without RP
filtered_data <- IOD_mat_melt[IOD_mat_melt$Pathways != "RP genes", ]
levene_result_noRP <- leveneTest(IOD ~ Pathways, data = filtered_data)

# Capture the result as a character string
result_text <- paste(capture.output(levene_result_noRP), collapse = "\n")

# Write the result to a text file with a message
write(paste("\n\nLevene Test Result UEG without RP:\n", result_text), 
      file = paste0(Tissue_name,"_","levene_test_result.txt"), 
      append = TRUE)

# pairwise comparisions to see if variance of RP is different from other genes
pairwise_results <- list()
groups <- unique(IOD_mat_melt$Pathways)

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    subset_data <- subset(IOD_mat_melt, Pathways %in% c(groups[i], groups[j]))
    pairwise_results[[paste(groups[i], groups[j], sep = "-")]] <- leveneTest(IOD ~ Pathways, data = subset_data)
  }
}

# Apply Bonferroni correction
adjusted_p_values <- p.adjust(sapply(pairwise_results, function(x) x$`Pr(>F)`[1]), method = "bonferroni")
adjusted_p_values <- as.data.frame(adjusted_p_values)

filename <- paste(Tissue_name, "Levenes_test_RP_Vs_UEG_pairwise",".csv", sep = "")
write.csv(adjusted_p_values,filename)


#RP Vs Random
#analysis 4 - comparision with random genes
identslist <- levels(seurat)
identslist <- as.list(identslist)
total_genes <- rownames(seurat)
genewithoutrp <- as.list(setdiff(total_genes, rpgenes))
# Number of sets to generate
num_sets <- 20

# List to store the sets of random genes
random_gene_sets <- list()

# Loop to generate random gene sets
for (i in 1:num_sets) {
  # Get a random sample of 100 genes
  random_genes <- as.list(sample(genewithoutrp, 76, replace = FALSE))
  
  #change it to list
  combined_random_genes <- unlist(random_genes, use.names = FALSE)
  
  # Store the random sample in the list
  random_gene_sets[[i]] <- combined_random_genes
}

write.csv(random_gene_sets,paste(Tissue_name,"Random_gene_set_sperm.csv"))


mat <- seurat[["SCT"]]$data[rpgenes, ] %>% as.matrix()
mean <- apply(mat, 1, mean) %>% as.numeric()
variance <- apply(mat, 1, var) %>% as.numeric()
final_variance_matrix_random <- as.matrix(variance)
rownames(final_variance_matrix_random) <- rpgenes
var_mat_melt <- melt(final_variance_matrix_UEG)
var_mat_melt$Var2 <- "RP genes"
IOD <- variance/mean
IOD_matrix <- as.matrix(IOD)
finalIOD_matrix_random <- IOD_matrix
rownames(finalIOD_matrix_random) <- rpgenes
IOD_mat_melt <- melt(finalIOD_matrix_random)
IOD_mat_melt$Var2 <- "RP genes"

for (j in 1:length(random_gene_sets)) {
  
  genelist <- as.list(random_gene_sets[j])
  genelist <- unlist(genelist)
  mat <- seurat[["SCT"]]$data[genelist, ] %>% as.matrix()
  mean <- apply(mat, 1, mean) %>% as.numeric()
  variance <- apply(mat, 1, var) %>% as.numeric()
  var_matrix <- as.matrix(variance)
  rownames(var_matrix) <- genelist
  var_mat_melt_j <- melt(var_matrix)
  var_mat_melt_j$Var2 <- j
  
  IOD <- variance/mean
  IOD_matrix <- as.matrix(IOD)
  rownames(IOD_matrix) <- genelist
  IOD_mat_melt_j <- melt(IOD_matrix)
  IOD_mat_melt_j$Var2 <- j
  
  var_mat_melt <- rbind(var_mat_melt,var_mat_melt_j)
  IOD_mat_melt <- rbind(IOD_mat_melt,IOD_mat_melt_j)
}

colnames(IOD_mat_melt) <- c("Genes","Gene_sets","IOD")
colnames(var_mat_melt) <- c("Genes","Gene_sets","IOD")

IOD_mat_melt$Gene_sets <- factor(IOD_mat_melt$Gene_sets, levels = c("RP genes","1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11","12","13","14","15","16","17","18","19","20"))

plot_random <- ggplot(IOD_mat_melt, 
                      aes(x = Gene_sets, y = IOD, fill = Gene_sets)) + 
  geom_violin(size = 0.1) + # Adjusts the outline thickness
  ylab("Index of Dispersion") +
  xlab("Genes") +
  ggtitle("Comparison of Index of dispersion of RP genes VS Random gene sets") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) +
  scale_fill_viridis_d(option = "G")  # Discrete color mapping

plot_random

file_path <- paste(Tissue_name, "RP Vs Random",".png", sep = "")
# Open a PNG device with the desired resolution
png(file = file_path, width = 10, height = 8, units = 'in', res = 600)
# Print the plot
plot(plot_random)
# Close the PNG device
dev.off()

var_mat_melt$Gene_sets <- factor(var_mat_melt$Gene_sets, levels = c("RP genes","1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11","12","13","14","15","16","17","18","19","20"))

#levene's test to test overall difference
levene_result_random <- leveneTest(IOD ~ Gene_sets, data = IOD_mat_melt)

# Capture the result as a character string
result_text <- paste(capture.output(levene_result_random), collapse = "\n")

# Write the result to a text file with a message
write(paste("\n\nLevene Test Result RP Vs Random:\n", result_text), 
      file = paste0(Tissue_name,"_","levene_test_result.txt"), 
      append = TRUE)

#Filtering the data to remove RPs and levene's test to see comparisions without RP
filtered_data <- IOD_mat_melt[IOD_mat_melt$Gene_sets != "RP genes", ]
levene_random_nonRP <- leveneTest(IOD ~ Gene_sets, data = filtered_data)

# Capture the result as a character string
result_text <- paste(capture.output(levene_random_nonRP), collapse = "\n")

# Write the result to a text file with a message
write(paste("\n\nLevene Test Result Random without RP:\n", result_text), 
      file = paste0(Tissue_name,"_","levene_test_result.txt"), 
      append = TRUE)

# pairwise comparisions to see if variance of RP is different from other genes
pairwise_results <- list()
groups <- unique(IOD_mat_melt$Gene_sets)

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    subset_data <- subset(IOD_mat_melt, Gene_sets %in% c(groups[i], groups[j]))
    pairwise_results[[paste(groups[i], groups[j], sep = "-")]] <- leveneTest(IOD ~ Gene_sets, data = subset_data)
  }
}

# Apply Bonferroni correction
adjusted_p_values <- p.adjust(sapply(pairwise_results, function(x) x$`Pr(>F)`[1]), method = "bonferroni")
adjusted_p_values <- as.data.frame(adjusted_p_values)

filename <- paste(Tissue_name, "Levenes_test_RP_Vs_random_pairwise",".csv", sep = "")
write.csv(adjusted_p_values,filename)


