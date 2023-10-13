
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
m6A_Pan <- read.csv("./support_files/ssGSEA.csv" , header = T, row.names = 1) ###
gene_set<-read.csv("./support_files/mmc3.csv")[, 1:2] ##
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])

####GSEA####
m6A_Pan_gsva_matrix<- gsva(as.matrix(m6A_Pan), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

Sample_Type <- read.csv("./support_files/SampleType.csv",header = T,row.names = 1)  ###

m6A_Pan_gsva_matrix<- t(scale(t(m6A_Pan_gsva_matrix)))
##########
library(pheatmap)

bk = unique(c(seq(-1.5,1.5, length=100))) ####

pdf(file = "./output_files/Pan_m6a_ssgsea.pdf",height = 6,width = 12)

pheatmap(m6A_Pan_gsva_matrix,
         show_colnames = F,scale = "row",border = F,
         cluster_rows = T,cluster_cols = F,
         breaks=bk,annotation_col = Sample_Type,
         treeheight_row = 20,cellheight=20,cellwidth=10,
         fontsize=10,gaps_col = c(10,14,17,19,23,31,37,42),width = 12)
dev.off()

