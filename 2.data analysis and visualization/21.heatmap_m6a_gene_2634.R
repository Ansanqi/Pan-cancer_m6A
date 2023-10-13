
library(pheatmap)
library(dplyr)

heatmap_m6a2607 <- read.csv("./support_files/m6A2607forheatmap.csv",header = T,row.names = 1)  ###
heatmap_m6a2607 <- heatmap_m6a2607[,c(1:84)]

Sample_Type <- read.csv("./support_files/SampleType84.csv",header = T)  ###
annotation_col = data.frame(Group = Sample_Type$Group,Type = Sample_Type$Type)#
rownames(annotation_col) = Sample_Type$Sample

heatmap_gene2607 <- read.csv("./support_files/gene2607forheatmap.csv",header = T,row.names = 1)  ###

gene_or <-as.matrix( row.names(heatmap_m6a2607)) ;colnames(gene_or)[1]<- "Gene"

heatmap_gene2607_or <- cbind(gene_or, heatmap_gene2607[match(gene_or[,1],row.names(heatmap_gene2607)),1])
heatmap_gene2607_or <- as.data.frame(heatmap_gene2607_or);b <- as.numeric(heatmap_gene2607_or[,2] );heatmap_gene2607_or <- cbind(heatmap_gene2607_or,b);heatmap_gene2607_or <- heatmap_gene2607_or[,-2];rm(b)
colnames(heatmap_gene2607_or)[1] <-"Gene";colnames(heatmap_gene2607_or)[2] <-"Colon_cancer_1";rm(gene_or)

for (i in 2:ncol(heatmap_gene2607)) {
  heatmap_gene2607_or <- cbind(heatmap_gene2607_or, heatmap_gene2607[match(heatmap_gene2607_or[,1],row.names(heatmap_gene2607)),i])
  colnames( heatmap_gene2607_or )[i+1]<-  colnames( heatmap_gene2607 )[i]
}   ###

row.names(heatmap_gene2607_or)<-heatmap_gene2607_or[,1];heatmap_gene2607_or <-heatmap_gene2607_or[,-1]

heatmap_gene2607_or <- log2(heatmap_gene2607_or+1)

pdf(file = "./output_files/heatmap_2607m6a_2607Gene.pdf")

pheatmap(heatmap_m6a2607,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = TRUE,cluster_cols = FALSE,
         show_rownames = FALSE,color=colorRampPalette(c("darkgreen","white","darkorange"))(1000))   #

pheatmap(heatmap_gene2607_or,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = FALSE,cluster_cols = FALSE,
         show_rownames = FALSE,color=colorRampPalette(c("darkgreen","white","darkorange"))(1000))   #

dev.off()


