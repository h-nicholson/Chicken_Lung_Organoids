####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#GOAL: Make data a Seurat Object

{
setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("hdf5r")

snRNAseq_invivo_h5 <- Read10X_h5("raw_data/snRNA-1/filtered_feature_bc_matrix.h5")
snRNAseq_invivo_matrix<-as.matrix(snRNAseq_invivo_h5)
snRNAseq_invivo_seurat <- CreateSeuratObject(counts = snRNAseq_invivo_h5)
snRNAseq_invivo_seurat[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_invivo_seurat, pattern = "^MT-")
snRNAseq_invivo_seurat <- subset(snRNAseq_invivo_seurat, subset = nFeature_RNA >= 200 & percent.mt <= 12)
snRNAseq_invivo_seurat <- NormalizeData(snRNAseq_invivo_seurat)

hist(log10(snRNAseq_invivo_matrix+1),breaks=100,ylim=c(0,10^5))




snRNAseq_org_h5 <- Read10X_h5("raw_data/snRNA-2/filtered_feature_bc_matrix.h5")
snRNAseq_org_matrix<-as.matrix(snRNAseq_org_h5)
snRNAseq_org_seurat <- CreateSeuratObject(counts = snRNAseq_org_h5)
snRNAseq_org_seurat[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_org_seurat, pattern = "^MT-")
snRNAseq_org_seurat <- subset(snRNAseq_org_seurat, subset = nFeature_RNA >= 200 & percent.mt <= 12)
snRNAseq_org_seurat <- NormalizeData(snRNAseq_org_seurat)


hist(log10(snRNAseq_org_matrix+1),breaks=100,ylim=c(0,10^5))



save(list=c("snRNAseq_org_seurat","snRNAseq_invivo_seurat"),file="snRNAseq_Seurat.RData")

}


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ##     ##     ########   ##########   #########   ########     ########    ###     ##      #########                                                      
##     ##          ##         ##     ##    ##              ##       ##          ##     ##       ##       ####    ##     ## 
##    ##           ##         ##     ##    ##              ##       ##          ##     ##       ##       ## ##   ##    ##                 
##    ##           ##         ##     ##     #######        ##       #######     ########        ##       ##  ##  ##    ##    ####
##    ##           ##         ##     ##           ##       ##       ##          ##    ##        ##       ##   ## ##    ##       ##
##     ##          ##         ##     ##           ##       ##       ##          ##     ##       ##       ##    ####     ##      ##       
##      ########   ########    #######     ########        ##       #########   ##      ##   ########    ##     ###      ########     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   #######   ##         ##               ##########  ##      ##  #######    ########    ########                                                                                         
##     ##          ##        ##         ##                   ##       ##    ##   ##    ##   ##         ##                 
##    ##           ##        ##         ##                   ##        ##  ##    ##    ##   ##         ##                                 
##    ##           ######    ##         ##        #####      ##         ####     #######    #######     #######                                
##    ##           ##        ##         ##                   ##          ##      ##         ##                ##      
##     ##          ##        ##         ##                   ##          ##      ##         ##                ##            
##      ########   #######   ########   ########             ##          ##      ##         ########    #######                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################



setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-11-22 - snRNAseq-chicken/")

library("Seurat")
library("ggplot2")

load("snRNAseq_Seurat.RData")


########################################################################################################################################
########################################################################################################################################
#1 IN VIVO DATA PREP::
########################################################################################################################################
########################################################################################################################################


{
seurat_object<-snRNAseq_invivo_seurat
########################################################################################################################################
#1 CLUSTER scRNAseq Data (unsupervised)
########################################################################################################################################


# Identify variable features across cells
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_object <- ScaleData(seurat_object)

# Perform PCA for dimensionality reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Find neighbors and clusters
seurat_object <- FindNeighbors(seurat_object)  # Adjust dimensions as needed
seurat_object <- FindClusters(seurat_object, resolution = 0.15)  # Adjust resolution as needed

# Run UMAP for visualization
seurat_object <- RunUMAP(seurat_object,dims=1:20)

# Plot the clusters using UMAP
DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("in vivo chicken snRNAseq")


########################################################################################################################################
#2 HEATMAP:  top 25 genes per cluster
########################################################################################################################################




library(dplyr)

# Ensure clustering has been performed and clusters are stored in 'seurat_clusters'

# Identify the top 10 variable genes for each cluster
top_genes <- seurat_object %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

# Extract the gene names
top_gene_list <- unique(top_genes$gene)[-grep("ENSGAL",unique(top_genes$gene))]

# Scale the data for heatmap visualization
seurat_object <- ScaleData(seurat_object, features = top_gene_list)

# Generate a heatmap
DoHeatmap(seurat_object, features = top_gene_list, group.by = "seurat_clusters") +
  ggtitle("Top 10 Variable Genes per Cluster")+
  theme(axis.text.y = element_text(size = 3)) 

########################################################################################################################################
#03 celltypes for each cluster
########################################################################################################################################


cluster="10"
top_genes$gene[which(top_genes$cluster==cluster)][-grep("ENSGAL",top_genes$gene[which(top_genes$cluster==cluster)])]
#CLUSTER 0: Endothelial cells (VWF, CDH5, EPAS1, EFNB2, )
#CLUSTER 1: Fibroblasts
#CLUSTER 2: ????                sensory neurons???
#CLUSTER 3: Epithelial         -> respiratory tract, basal, tight junctions
#CLUSTER 4: Fibroblasts (PVLs)
#CLUSTER 5: Tcells 
#CLUSTER 6:  Proliferation      (epithelial, hematopoetic progenetors)
#CLUSTER 7: RBC's                 (HBAA, HBAD, RAP1GAP2)
#CLUSTER 8: Smooth muscle               (SYNPO2, STAC, MYOCD, COL19A1, MYH11, CNN1, NEXN, RYR3)
#CLUSTER 9: ciliated cells (RSPH14, )
#CLUSTER 10: alvelolar cells or secretory (ELF5, AQP3, )



########################################################################################################################################
#4 CHECK UNSUPERVISED ASSIGNMENT BASED ON CANONICAL MARKERS: (check if roughly consistent)
########################################################################################################################################

seurat_4cells<-seurat_object


# Single-Cell Transcriptome Atlas of Newcastle Disease Virus in Chickens Both In Vitro and In Vivo
markers <- list(fibroblast = c("COLLA2","MGP","FN1","PDGFRa","FHl1","LOX","LIMHC1","MFAP1","MFAP2","MFAP5","MFAP3L","TCF21","TCF3","TCF12","TCF25","TCF7L2","HSD11B2"),
                epithelial = c("HOPX","SFTPA1","SFTPA2","TPPP3","KRT18","KRT5","KRT7","KRT12","KRT13","KRT14","KRT15","KRT18","KRT19","CLIC5","CTSH","CLDN18","Cbr1","Cbr4"),
                endothelial = c("GJA5","BMX","CA4","PROX1","CAV1","CAV2","CAV3","CDH1","CDH2","CDH3","CDH4","CDH5","CDH6","CD93","LYVE1","PTPRB","TEK","KDR","NRP1","FIBIN","CLEC3A","CLEC3B","CLEC16A","CLEC19A","EDNRB","TBX1","TBX2","TBX3","TBX4","TBX5","TBX6","TIMP3"),
                tcell = c("CD3D","CD3E","CD4","CD8A","Thy1","LEF1","CD28","ICOS","CD8A","CCR6","LEF1","IL7R"),
                myeloid = c("ITGB2","CST3","CST7","CSF1R","SIRPA","NAAA","CSF3R","S100A9","EGR1","CD14","TREM2","CFD"))

#Markers for all cell-types
markers <- list(# Endothelial cells
  endothelial = c("VWF", "PECAM1", "CDH5", "VEGFR2", "NOS3", "ICAM1", "EFNB2", "FLT1", "CLDN5", "SELE", 
                         "ANGPT2", "TIE1", "KDR", "PLVAP", "EHD2", "CD34", "ESAM", "JAG1", "CCL21", "MMP2"),
  fibroblast = c("COL1A1", "COL3A1", "FBLN1", "FBLN5", "PDGFRA", "ACTA2", "THY1", "LUM", "DCN", "SPARC", 
                        "SFRP4", "LOX", "MFAP4", "POSTN", "TAGLN", "PDPN", "P4HA1", "ELN", "FN1", "ITGB1"),
  basal = c("TP63", "KRT5", "KRT14", "KRT15", "ITGA6", "ITGB4", "CD44", "NGFR", "BCL11A", "LGR5", 
                              "EGFR", "S100A2", "SOX2", "SERPINB5", "CXCL14", "MMP9", "PKP1", "COL17A1", "DSP", "GJB2"),
  tcell = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4", "IL2RA", "FOXP3", "CCR7", "CD28", "IKZF1", 
                   "CTLA4", "GZMB", "IFNG", "IL7R", "PDCD1", "CD69", "CD44", "NFATC2", "STAT5A", "CXCR3"),
  prolif = c("MKI67", "PCNA", "TOP2A", "AURKA", "BIRC5", "CCNB1", "CDK1", "UBE2C", "CENPF", "CENPA", 
                           "NDC80", "BUB1", "MAD2L1", "PLK1", "CDC20", "E2F1", "RRM2", "TK1", "CHEK1", "MCM2"),
  RBC = c("HBA1", "HBA2", "HBB", "AHSP", "EPOR", "KLF1", "GATA1", "BCL11A", "ANK1", "SLC4A1", 
                            "ALAS2", "CA1", "CA2", "ADD2", "RHCE", "RHAG", "GYPA", "GYPC", "BPGM", "TRIM58"),
  muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1", "LMOD1", "CAV1", "TPM2", "TPM1", "CALD1", "SMTN", 
                           "MYLK", "FHL1", "VIM", "DES", "MYOCD", "PLN", "PDGFRB", "NOTCH3", "ELN", "KCNMB1"),
  ciliated = c("FOXJ1", "DNAH5", "DNAH9", "DNAI1", "RSPH1", "TEKT1", "CFAP43", "CFAP69", "CCDC39", "CCDC40", 
                                 "TPPP3", "PIFO", "SNTN", "SPAG1", "ARMC4", "DNAAF1", "LRRC6", "IFT81", "SPAG17", "ODF2"),
  secretory = c("MUC1", "MUC4", "MUC5B", "MUC16", "PIGR", "PRLR", "ELF5", "SPDEF", "SCGB1A1", "WFDC2", 
                                  "AQP5", "LTF", "LALBA", "KRT7", "KRT18", "AGR2", "CLCA1", "EPHB2", "SERPINB3", "LCN2")
  
)



# Score each cell type based on the defined markers
for (cell_type in names(markers)) {
  seurat_4cells <- AddModuleScore(
    object = seurat_4cells, 
    features = list(markers[[cell_type]]), 
    name = paste0(cell_type, "_Score")
  )
}

# Assign cell type identities based on the highest module score for each cell
cell_type_scores <- as.data.frame(seurat_4cells@meta.data[, grep("_Score", colnames(seurat_4cells@meta.data))])
cell_type_assignment <- colnames(cell_type_scores)[max.col(cell_type_scores, ties.method = "first")]

# Add the cell type assignment to the Seurat object metadata
seurat_4cells$cell_type <- gsub("_Score", "", cell_type_assignment)


colors<-c("#B22222", "#1F78B4", "#33A02C", "darkorange","purple","black","grey","darkcyan","yellow4")
# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_4cells, reduction = "umap", group.by = "cell_type") +
  ggtitle("Mouse Retina Cell Type Assignment based on Marker Genes")+
  scale_color_manual(values=colors)


########################################################################################################################################
#5 Assign FINAL CELL TYPE: by cluster
########################################################################################################################################

tmp_cell_ids<-seurat_object$seurat_clusters

tmp_cell_ids<-gsub("0","Endothelial",tmp_cell_ids)
tmp_cell_ids<-gsub("1","Fibroblasts",tmp_cell_ids)
tmp_cell_ids<-gsub("2","Epithelial",tmp_cell_ids)
tmp_cell_ids<-gsub("3","Epithelial",tmp_cell_ids)
tmp_cell_ids<-gsub("4","Fibroblasts",tmp_cell_ids)
tmp_cell_ids<-gsub("5","T-cells",tmp_cell_ids)
tmp_cell_ids<-gsub("6","Prolifer?",tmp_cell_ids)
tmp_cell_ids<-gsub("7","RBC",tmp_cell_ids)
tmp_cell_ids<-gsub("8","Muscle",tmp_cell_ids)
tmp_cell_ids<-gsub("9","Ciliated",tmp_cell_ids)
tmp_cell_ids<-gsub("FibroblastsEndothelial","Aveloar/Secr?",tmp_cell_ids)


seurat_object$cell_type <- tmp_cell_ids


colors<-c("darkcyan","blue","darkred","darkblue","grey40","grey","darkorange","red","darkgreen")
# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_object, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo chicken lung")+
  scale_color_manual(values=colors)


snRNAseq_invivo_cellTypes<-seurat_object

}



########################################################################################################################################
########################################################################################################################################
#2 ORGANOID DATA PREP::
########################################################################################################################################
########################################################################################################################################

{
seurat_object<-snRNAseq_org_seurat
########################################################################################################################################
#1 CLUSTER scRNAseq Data (unsupervised)
########################################################################################################################################


# Identify variable features across cells
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_object <- ScaleData(seurat_object)

# Perform PCA for dimensionality reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Find neighbors and clusters
seurat_object <- FindNeighbors(seurat_object)  # Adjust dimensions as needed
seurat_object <- FindClusters(seurat_object, resolution = 0.15)  # Adjust resolution as needed

# Run UMAP for visualization
seurat_object <- RunUMAP(seurat_object,dims=1:20)

# Plot the clusters using UMAP
DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("in vivo chicken snRNAseq")


########################################################################################################################################
#2 HEATMAP:  top 25 genes per cluster
########################################################################################################################################




library(dplyr)

# Ensure clustering has been performed and clusters are stored in 'seurat_clusters'

# Identify the top 10 variable genes for each cluster
top_genes <- seurat_object %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

# Extract the gene names
top_gene_list <- unique(top_genes$gene)[-grep("ENSGAL",unique(top_genes$gene))]

# Scale the data for heatmap visualization
seurat_object <- ScaleData(seurat_object, features = top_gene_list)

# Generate a heatmap
DoHeatmap(seurat_object, features = top_gene_list, group.by = "seurat_clusters") +
  ggtitle("Top 10 Variable Genes per Cluster")+
  theme(axis.text.y = element_text(size = 3)) 

########################################################################################################################################
#03 celltypes for each cluster
########################################################################################################################################

cluster="2"
top_genes$gene[which(top_genes$cluster==cluster)][-grep("ENSGAL",top_genes$gene[which(top_genes$cluster==cluster)])]
#CLUSTER 0: ?????????           neurons?
#CLUSTER 1: ????????            DEAD??  ribosomes, mitochondria?   ->. progenitor cells, stem cells, macrophages/DC's

#CLUSTER 2: proliferative
#CLUSTER 3: fibroblasts


FeaturePlot(seurat_object, features = "percent.mt", cols = c("lightblue", "darkred"))


seurat_object@meta.data$read_depth<-log10(seurat_object@meta.data$nFeature_RNA+1)

FeaturePlot(seurat_object, features = "read_depth", cols = c("lightblue", "darkred"))


DimPlot(seurat_object, group.by = "nCount_RNA") +
  scale_color_gradient(low = "lightblue", high = "darkred")


########################################################################################################################################
#4 CHECK UNSUPERVISED ASSIGNMENT BASED ON CANONICAL MARKERS: (check if roughly consistent)
########################################################################################################################################

seurat_4cells<-seurat_object
markers <- list(
  # Marker genes for Epithelial Cells
  epithelial = c("EPCAM", "KRT18", "KRT19", "CDH1", "MUC1", "OCLN", "CLDN3", "CLDN4", "CLDN7", "TJP1", "VIM", "SFTPC", "SFTPA1", "SFTPA2", "KRT5", "KRT14", "TP63", "FOXJ1", "MUC5AC", "MUC5B"),
  
  # Marker genes for Fibroblasts
  fibroblast = c("PDGFRA", "COL1A1", "ACTA2", "FAP", "VIM", "COL3A1", "FN1", "THY1", "SERPINF1", "MMP2", "PDGFRB", "POSTN", "S100A4", "ITGB1", "DCN", "SPARC", "COL5A1", "ITGA11", "LUM", "FMOD")
)



seurat_4cells<-seurat_object
markers <- list(
  # Marker genes for Alveolar Cells (Type I and Type II combined)
  alveolar_genes = c("AGER", "PDPN", "HOPX", "CAV1", "CLDN18", "SPP1", "CAV2", "S1PR1", "RTKN2", "AQP5", "C1QTNF3", "T1alpha", "RAGE", "AQP3", "CYP2F2", "CYP2A1", "GPRC5A", "FABP5", "RGS5", "PECAM1"),
  
  # Marker genes for Basal Epithelial Cells
  basal_genes = c("KRT5", "TP63", "KRT14", "NGFR", "CD44", "ITGA6", "JAG1", "CDH3", "NCAM1", "PDPK1", "CDH2", "KRT17", "CDH1", "TGFB1", "GATA6", "SOX2", "MUC5B", "NOTCH3", "FGFR2", "HES1"),
  
  # Marker genes for Ciliated Cells
  ciliated_genes = c("FOXJ1", "DNAI1", "TEKT1", "HYDIN", "DNAH5", "SPAG6", "RSPH9", "MCIDAS", "GAS8", "PIFO", "CCDC39", "CCDC40", "DRC1", "DNALI1", "CFAP54", "ARMC4", "SPAG16", "CFAP300", "HYDIN2", "C2CD3"),
  
  # Marker genes for Fibroblasts
  fibroblast_genes = c("PDGFRA", "COL1A1", "ACTA2", "FAP", "VIM", "COL3A1", "FN1", "THY1", "SERPINF1", "MMP2", "PDGFRB", "POSTN", "S100A4", "ITGB1", "DCN", "SPARC", "COL5A1", "ITGA11", "LUM", "FMOD")
)





seurat_4cells<-seurat_object
markers <- list(# Endothelial cells
  endothelial = c("VWF", "PECAM1", "CDH5", "VEGFR2", "NOS3", "ICAM1", "EFNB2", "FLT1", "CLDN5", "SELE", 
                  "ANGPT2", "TIE1", "KDR", "PLVAP", "EHD2", "CD34", "ESAM", "JAG1", "CCL21", "MMP2"),
  fibroblast = c("COL1A1", "COL3A1", "FBLN1", "FBLN5", "PDGFRA", "ACTA2", "THY1", "LUM", "DCN", "SPARC", 
                 "SFRP4", "LOX", "MFAP4", "POSTN", "TAGLN", "PDPN", "P4HA1", "ELN", "FN1", "ITGB1"),
  basal = c("TP63", "KRT5", "KRT14", "KRT15", "ITGA6", "ITGB4", "CD44", "NGFR", "BCL11A", "LGR5", 
            "EGFR", "S100A2", "SOX2", "SERPINB5", "CXCL14", "MMP9", "PKP1", "COL17A1", "DSP", "GJB2"),
  tcell = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4", "IL2RA", "FOXP3", "CCR7", "CD28", "IKZF1", 
            "CTLA4", "GZMB", "IFNG", "IL7R", "PDCD1", "CD69", "CD44", "NFATC2", "STAT5A", "CXCR3"),
  prolif = c("MKI67", "PCNA", "TOP2A", "AURKA", "BIRC5", "CCNB1", "CDK1", "UBE2C", "CENPF", "CENPA", 
             "NDC80", "BUB1", "MAD2L1", "PLK1", "CDC20", "E2F1", "RRM2", "TK1", "CHEK1", "MCM2"),
  RBC = c("HBA1", "HBA2", "HBB", "AHSP", "EPOR", "KLF1", "GATA1", "BCL11A", "ANK1", "SLC4A1", 
          "ALAS2", "CA1", "CA2", "ADD2", "RHCE", "RHAG", "GYPA", "GYPC", "BPGM", "TRIM58"),
  muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1", "LMOD1", "CAV1", "TPM2", "TPM1", "CALD1", "SMTN", 
             "MYLK", "FHL1", "VIM", "DES", "MYOCD", "PLN", "PDGFRB", "NOTCH3", "ELN", "KCNMB1"),
  ciliated = c("FOXJ1", "DNAH5", "DNAH9", "DNAI1", "RSPH1", "TEKT1", "CFAP43", "CFAP69", "CCDC39", "CCDC40", 
               "TPPP3", "PIFO", "SNTN", "SPAG1", "ARMC4", "DNAAF1", "LRRC6", "IFT81", "SPAG17", "ODF2"),
  secretory = c("MUC1", "MUC4", "MUC5B", "MUC16", "PIGR", "PRLR", "ELF5", "SPDEF", "SCGB1A1", "WFDC2", 
                "AQP5", "LTF", "LALBA", "KRT7", "KRT18", "AGR2", "CLCA1", "EPHB2", "SERPINB3", "LCN2")
  
)



# Score each cell type based on the defined markers
for (cell_type in names(markers)) {
  seurat_4cells <- AddModuleScore(
    object = seurat_4cells, 
    features = list(markers[[cell_type]]), 
    name = paste0(cell_type, "_Score")
  )
}

# Assign cell type identities based on the highest module score for each cell
cell_type_scores <- as.data.frame(seurat_4cells@meta.data[, grep("_Score", colnames(seurat_4cells@meta.data))])
cell_type_assignment <- colnames(cell_type_scores)[max.col(cell_type_scores, ties.method = "first")]

# Add the cell type assignment to the Seurat object metadata
seurat_4cells$cell_type <- gsub("_Score", "", cell_type_assignment)


colors<-c("#B22222", "#1F78B4", "#33A02C", "darkorange","purple","black","grey","darkcyan","yellow4")
# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_4cells, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Types:  organoid chicken lung")+
  scale_color_manual(values=colors)


########################################################################################################################################
#5 Assign FINAL CELL TYPE: by cluster
########################################################################################################################################

tmp_cell_ids<-seurat_object$seurat_clusters

tmp_cell_ids<-gsub("0","Epithelial",tmp_cell_ids)
tmp_cell_ids<-gsub("1","Low quality?",tmp_cell_ids)
tmp_cell_ids<-gsub("2","Epithelial (cycling)",tmp_cell_ids)
tmp_cell_ids<-gsub("3","Fibroblasts",tmp_cell_ids)


seurat_object$cell_type <- tmp_cell_ids


colors<-c("darkblue","blue","grey40","grey90")
# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_object, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo chicken lung")+
  scale_color_manual(values=colors)


snRNAseq_organd_cellTypes<-seurat_object


save(list=c("snRNAseq_organd_cellTypes","snRNAseq_invivo_cellTypes"),file="snRNAseq_cellTypes.RData")


}



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ########        ##       ###     ##                                  
##     ##          ##         ##             ####      ####    ##                             
##    ##           ##         ##            ##  ##     ## ##   ##                                             
##    ##           ##         #######      ##    ##    ##  ##  ##                                  
##    ##           ##         ##          ##########   ##   ## ##                                   
##     ##          ##         ##          ##      ##   ##    ####                                   
##      ########   ########   ########    ##      ##   ##     ###                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#DATE:  1-3-2024
#GOAL:  they want boiler plate seurat analysis




setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-11-22 - snRNAseq-chicken/")

library("Seurat")
library("ggplot2")

load("snRNAseq_cellTypes_CLEAN.RData")

####################################################################################################################################################################################################################################################################
#CLEAN DATA:
#################################################################################################################################################################################################################################################################### 

{
load("snRNAseq_cellTypes.RData")
################################################################
#CHECK IN VIVO
################################################################
snRNAseq_invivo_cellTypes@meta.data$read_depth<-log10(snRNAseq_invivo_cellTypes@meta.data$nFeature_RNA+1)
FeaturePlot(snRNAseq_invivo_cellTypes, features = "read_depth", cols = c("lightblue", "darkred"))



################################################################
#CLEAN IN VITRO:
################################################################

#CHECK IN VITRO"
FeaturePlot(snRNAseq_organd_cellTypes, features = "read_depth", cols = c("lightblue", "darkred"))

DimPlot(snRNAseq_organd_cellTypes, reduction = "umap", group.by = "cell_type")
seurat1<-subset(snRNAseq_organd_cellTypes, subset = cell_type != "Low quality?")
seurat1 <- subset(seurat1, subset = nFeature_RNA >= 2500 & percent.mt <= 12)

DimPlot(seurat_object, reduction = "umap", group.by = "cell_type")
FeaturePlot(seurat1, features = "read_depth", cols = c("lightblue", "darkred"))


seurat1 <- NormalizeData(seurat1)

seurat_object<-seurat1

# Identify variable features across cells
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_object <- ScaleData(seurat_object)

# Perform PCA for dimensionality reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Find neighbors and clusters
seurat_object <- FindNeighbors(seurat_object)  # Adjust dimensions as needed
seurat_object <- FindClusters(seurat_object, resolution = 0.15)  # Adjust resolution as needed

# Run UMAP for visualization
seurat_object <- RunUMAP(seurat_object,dims=1:20)

# Plot the clusters using UMAP
DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("in vivo chicken snRNAseq")


################################
#2 HEATMAP:  top 25 genes per cluster
################################




library(dplyr)

# Ensure clustering has been performed and clusters are stored in 'seurat_clusters'

# Identify the top 10 variable genes for each cluster
top_genes <- seurat_object %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

# Extract the gene names
top_gene_list <- unique(top_genes$gene)[-grep("ENSGAL",unique(top_genes$gene))]

# Scale the data for heatmap visualization
seurat_object <- ScaleData(seurat_object, features = top_gene_list)

# Generate a heatmap
DoHeatmap(seurat_object, features = top_gene_list, group.by = "seurat_clusters") +
  ggtitle("Top 10 Variable Genes per Cluster")+
  theme(axis.text.y = element_text(size = 3)) 

################################
#03 celltypes:  choose not to name as ambiguous still...
################################

cluster="3"
top_genes$gene[which(top_genes$cluster==cluster)][-grep("ENSGAL",top_genes$gene[which(top_genes$cluster==cluster)])]
#CLUSTER 0: neuroendocrine?/
#CLUSTER 1: Aveloar type 2??

#CLUSTER 2: proliferative
#CLUSTER 3: fibroblasts




snRNAseq_organd_cellTypes<-seurat_object


save(list=c("snRNAseq_organd_cellTypes","snRNAseq_invivo_cellTypes"),file="snRNAseq_cellTypes_CLEAN.RData")

library(qs)
qsave(snRNAseq_organd_cellTypes, "shiny_orgnoid.qs")
qsave(snRNAseq_invivo_cellTypes, "shiny_invivo.qs")

#seurat_merged <- qs::qread("FINAL_B816-organoids.qs")



}



####################################################################################################################################################################################################################################################################
#BOILER PLATE:
#################################################################################################################################################################################################################################################################### 


# Plot the clusters using UMAP
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("in vivo chicken snRNAseq")



colors<-c("darkcyan","blue","darkred","darkblue","grey40","grey","darkorange","red","darkgreen")
# Visualize the assigned cell types in the UMAP plot
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo chicken lung")+
  scale_color_manual(values=colors)

seurat_obj<-snRNAseq_organd_cellTypes


colors<-c("darkblue","blue","grey40")
# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid chicken lung")+
  scale_color_manual(values=colors)

############################################################
#Cell-Type-oriented seurat plots
############################################################

library(Seurat)
library(tidyverse)


# Set "cell_type" as the grouping identity
Idents(seurat_obj) <- seurat_obj@meta.data$cell_type

# Exclude genes containing "ENSGAL"
valid_genes <- rownames(seurat_obj)[!grepl("ENSGAL", rownames(seurat_obj))]

# Subset Seurat object to exclude these genes
seurat_obj <- subset(seurat_obj, features = valid_genes)

# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(top_genes_per_cell_type, "organoid2_TopGenes.csv", row.names = FALSE)

# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("organoid2_DotPlot.pdf", plot = dotplot)

# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
})

pdf("organoid2_FeaturePlots.pdf", width = 10, height = 10)
for (plot in feature_plots) {
  print(plot)
}
dev.off()

# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("organoid2_Heatmap.pdf", plot = heatmap)


#############################
# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
#############################
library(patchwork) # For combining plots into a grid


top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

# Create individual feature plots for each gene
feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(gene) +
    theme_minimal()
})

# Combine all feature plots into a 5x2 grid
combined_plot <- wrap_plots(feature_plots, ncol = 2)

# Save the combined plot as a single PDF
ggsave("organoid2_FeaturePlots.pdf.pdf", plot = combined_plot, width = 15, height = 8)





