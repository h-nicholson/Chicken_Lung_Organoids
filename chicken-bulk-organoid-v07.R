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
#DATE:   11-22-2022
#GOAL:   Basic Parsing of GEO Dataset 
#PLAN:
setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-11-22 - snRNAseq-chicken/04 Bulk RNAseq/")

{

load("Chicken_RNAseq.RData")

#####################################################
#1.1 HISTOGRAMS:  make sure data is log-normally distributed
#####################################################

hist(log10(Chicken_counts+1),breaks=50,ylim=c(0,40000),main="Gene Expression is ~ Log-Normal Distributed")

#####################################################
#1.2 BARPLOTS:  make sure each sample has good "read depth"
#####################################################


read_depth<-colSums(Chicken_counts)

barplot(sort(read_depth),las=2,cex.names = 0.6,main="Read Depth = Amount mRNA/sample",)
abline(h=10e6,lty=2)


#####################################################
#1.3 NORMALIZE COLUMNS (read depth) 
#####################################################


barplot(sort(colSums(Chicken_tpm)),las=2,cex.names = 0.6,main="Gene Expression normalized to Depth")



#####################################################
#2.1 NORMALIZE ROWS (z-score)
#####################################################

#FILTER OUT ALL ZERO ROWS:
zero_rows<-which(rowSums(Chicken_tpm)==0)
organoid_logRPM_zero<-Chicken_tpm[-zero_rows,]

Chicken_tpm<-organoid_logRPM_zero

#CONVERT TO Z-SCORE
allSamples_zscore<-t(scale(t(Chicken_tpm),center=T,scale=T))


#####################################################
#NORMALIZE TISSUE AND ORGANOID SEPARALY
#####################################################
organoid_zscore<-t(scale(t(Chicken_tpm[,grep("orgnd",colnames(Chicken_tpm))]),center=T,scale=T))
tissue_zscore<-t(scale(t(Chicken_tpm[,grep("issue",colnames(Chicken_tpm))]),center=T,scale=T))

separated_organ_tissue_zscore_raw<-cbind(organoid_zscore,tissue_zscore)

separated_organ_tissue_zscore<-separated_organ_tissue_zscore_raw[-which(is.na(rowSums(separated_organ_tissue_zscore_raw))==T),]

#####################################################
#2.2 SAVE as CSV (to open in excel) and as RData file (so easy to do future R analysis on)
#####################################################


write.csv(allSamples_zscore,file="allSamples_zscore.csv")

save(list=c("allSamples_zscore","Chicken_tpm","separated_organ_tissue_zscore"),file="Tissue-Organoids_zscores.RData")

}  
  


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##       ##        ##        ##      ##   ########    ##     #######                               
##      ####     ####    ##      ####       ##         ##    ##   ##           ##    ##                                               
##     ##  ##    ## ##   ##     ##  ##      ##          ##  ##    ##           ##    ##           
##    ##    ##   ##  ##  ##    ##    ##     ##           ####      #######     ##     #######                                           
##    ########   ##   ## ##   ##########    ##            ##             ##    ##           ##                                      
##    ##    ##   ##    ####   ##      ##    ##            ##             ##    ##           ##                                          
##    ##    ##   ##     ###   ##      ##    ########      ##       #######     ##     #######                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE:   11-22-2022
#GOAL:   Basic Analysis of RNAseq dataset
#PLAN:

#       1.  UNSUPERVISED ANALYSIS:   
#                           -> CORRELATION MATRIX:   are replicates correlated?   what groups are there?
#                           -> PRINCIPAL COMPONENT ANALYSIS:   what are the 2 major "dimensions" of varaibility in samples
#                           -> PATHWAY ANALYSIS:  transcriptional programs
#
#################################################################################################################################################################################################################################################################### 

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-11-22 - snRNAseq-chicken/04 Bulk RNAseq/")


#load("Chicken_paired.RData")
load("Tissue-Organoids_zscores.RData")

library("gplots")
library("RColorBrewer")

###############################################################################################################################################################
#PART 1:   ORGANOIDS VS TISSUE
###############################################################################################################################################################


#####################################################
#FIGURE 2A:  CORR-MATRIX.    Published analysis
#####################################################

organoid_rpm_corr<-cor(Chicken_tpm)

heatmap.2(organoid_rpm_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "rpm corr",key.xlab = "rpm corr",cexRow = 0.6,cexCol = 0.6)

#####################################################
#FIGURE 2A:   CORR-MATRIX.   Refined analysis
#####################################################

organoid_zscore_corr<-cor(allSamples_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.6,cexCol = 0.6)





#####################################################
#FIGURE 2B:   PCA   published analysis
#####################################################

table(sapply(strsplit(colnames(allSamples_zscore),"_"), function(x) x[2]))
#Bladder Endometrium      Kidney       Liver        Lung       ovary       Ovary    Pancreas 
#4       4                4             4           4           1           1           6 

#RUN PCA
organoid_pca<-prcomp(t(allSamples_zscore))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)


sample_names_vector[grep("tissue_lung1",sample_names_vector)]<-"darkblue"
sample_names_vector[grep("orgnd_lung1",sample_names_vector)]<-"blue"


sample_names_vector[grep("tissue_lung7",sample_names_vector)]<-"red4"
sample_names_vector[grep("orgnd_lung7",sample_names_vector)]<-"red"

sample_names_vector[grep("tissue_lung3",sample_names_vector)]<-"darkgreen"
sample_names_vector[grep("orgnd_lung3",sample_names_vector)]<-"chartreuse4"



plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16,
     xlab="PC1 (54%):  organoid vs tissue",ylab="PC2 (11%):  Lung 7 vs 3/1")
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)

plot(organoid_pca_data[,1],organoid_pca_data[,4],col=sample_names_vector,pch=16,
     xlab="PC1 (54%):  organoid vs tissue",ylab="PC4 (5%): organoid = tissue")
text(organoid_pca_data[,1],organoid_pca_data[,4],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)


#####################################################
#heatmaps
#####################################################

#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Chicken_tpm,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(allSamples_zscore[genes_in_order[1:5000],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.1,cexCol = 0.5,mar=c(5,10))




#####################################################
#SELECT GENE SETS
#####################################################

cellType_genes<-read.csv("04 Requests/chicken lung cell type genes.csv",as.is=T,header=F)$V1
viral_genes<-read.csv("04 Requests/chicken lung viral genes.csv",as.is=T,header=F)$V1

#####################
#Cell-Type Genes:
#####################
cell_genes<-rownames(allSamples_zscore)[which(rownames(allSamples_zscore) %in% cellType_genes)]

#all samples z-score
heatmap.2(allSamples_zscore[cell_genes,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.5,cexCol = 0.5,mar=c(5,10))


#w/n cell-type z-score
heatmap.2(separated_organ_tissue_zscore[cell_genes,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.5,cexCol = 0.5,mar=c(5,10))


#tpm
heatmap.2(log10(Chicken_tpm+1)[cell_genes,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "logTPM",
          cexRow = 0.5,cexCol = 0.5,mar=c(5,10))


#####################
#Viral Genes:
#####################
viral_genes<-rownames(allSamples_zscore)[which(rownames(allSamples_zscore) %in% viral_genes)]

heatmap.2(allSamples_zscore[viral_genes,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.4,cexCol = 0.5,mar=c(5,10))



#w/n cell-type z-score
heatmap.2(separated_organ_tissue_zscore[viral_genes,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score (intnl)",
          cexRow = 0.5,cexCol = 0.5,mar=c(5,10))


#tpm
heatmap.2(log10(Chicken_tpm+1)[viral_genes,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "logTPM",
          cexRow = 0.5,cexCol = 0.5,mar=c(5,10))



#####################################################
#Venn Diagram
#####################################################

library("eulerr")

organoid_avg<-rowMeans(Chicken_tpm[,grep("orgnd",colnames(allSamples_zscore))])
tissue_avg<-rowMeans(Chicken_tpm[,grep("tissue",colnames(allSamples_zscore))])

threshold=0.1
scRNAseq_venn_list<-list(organoid=names(which(organoid_avg > threshold)),
                         tissue=names(which(tissue_avg > threshold)))



plot(euler(scRNAseq_venn_list, shape = "ellipse"), quantities = TRUE,main="Comparison of Expressed Genes (tpm)")

plot(log10(organoid_avg+1),log10(tissue_avg+1),pch=16,cex=0.5,xlab="Organoid Expression (logTPM)",ylab="Tissue Expression (logTPM)")
abline(a=0,b=1)


#####################################################
#FIGURE 3:   PATHWAY ANALYSIS
#####################################################

#####################
#LOAD PACKAGES
#####################


library(viper)
library(msigdbr)



#####################
#2 MAKE HALLMARK-OBJECT FOR PATHWAY ANALYSIS
#####################

#EXTRACT MOUSE-HALLMARKS:
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets$gs_name<-gsub("HALLMARK_","",HALLMARK_gene_sets$gs_name)


#CONVERT MOUSE-HALLMARK DATA-FRAME TO A "REGULON"
msigdbr_list = split(x = HALLMARK_gene_sets$gene_symbol, f = HALLMARK_gene_sets$gs_name)
hallmark_regulon<-list()

for (tf in names(msigdbr_list)){
  sig_genes<-msigdbr_list[[tf]]
  
  hallmark_regulon[[tf]]<-list()
  tmp_tfmode<-rep(1,length(sig_genes))
  names(tmp_tfmode)<-sig_genes
  
  hallmark_regulon[[tf]]$'tfmode'<-tmp_tfmode
  hallmark_regulon[[tf]]$'likelihood'<-as.numeric(tmp_tfmode)
}

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(allSamples_zscore,hallmark_regulon,method="scale")



#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.5,mar=c(5,10))



#####################################################
#FIGURE 3:   PATHWAY BARPLOT:
#####################################################



organ_vs_tissue<-sort(rowMeans(organoid_Hallmarks[,grep("orgnd",colnames(organoid_Hallmarks))])-rowMeans(organoid_Hallmarks[,grep("issue",colnames(organoid_Hallmarks))]))

barplot(organ_vs_tissue[c(1:15,40:50)],horiz=T,las=2,cex.names = 0.4,col="black",mar=c())



###############################################################################################################################################################
#PART 2:   ORGANOIDS and TISSUE
###############################################################################################################################################################




#RUN PCA
organoid_pca<-prcomp(t(separated_organ_tissue_zscore))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)


sample_names_vector[grep("tissue_lung1",sample_names_vector)]<-"darkblue"
sample_names_vector[grep("orgnd_lung1",sample_names_vector)]<-"blue"


sample_names_vector[grep("tissue_lung7",sample_names_vector)]<-"red4"
sample_names_vector[grep("orgnd_lung7",sample_names_vector)]<-"red"

sample_names_vector[grep("tissue_lung3",sample_names_vector)]<-"darkgreen"
sample_names_vector[grep("orgnd_lung3",sample_names_vector)]<-"chartreuse4"



plot(organoid_pca_data[,1],organoid_pca_data[,3],col=sample_names_vector,pch=16,
     xlab="PC1 (22%)",ylab="PC3 (11%)")
text(organoid_pca_data[,1],organoid_pca_data[,3],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)






sep_Zscore_Hallmarks <- viper(separated_organ_tissue_zscore,hallmark_regulon,method="scale")




#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(sep_Zscore_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(sep_Zscore_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.5,mar=c(5,10))

