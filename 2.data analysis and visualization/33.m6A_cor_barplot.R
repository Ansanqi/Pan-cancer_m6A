library(ggpubr)
library(reshape2)

cor_RBP<- read.csv("./support_files/cor_RBP.csv",header = T,row.names = 1)
gene_type<- read.csv("./support_files/gene_RBP.csv",header = T)
cor_RBP <- abs(cor_RBP)

sum <- rowSums(cor_RBP[,c(1:ncol(cor_RBP))])
sum <- t(sum);sum <- t(sum)

row_mean = apply(cor_RBP,1,mean)
row_mean <- t(row_mean);row_mean <- t(row_mean)


work_df <- cbind(gene_type,row_mean)
work_df <- cbind(gene_type,sum)
work_df <- as.data.frame(work_df)

pdf(file = "./output_files/boxplot.pdf")
ggboxplot(data = work_df, x = 'Description', y = 'sum', color = 'Description') +stat_compare_means(method = 't.test', comparisons = list(c('classical', 'referenced'), c('classical', 'extra'), c('referenced', 'extra')))#
dev.off()


m6A_84sam<- read.csv("./support_files/Allsample_m6A_new.csv",header = T,row.names = 1)
ens <- read.table("./support_files/Hg38ENSG.txt",header = T,row.names = 1)

setwd("./support_files/gene_abund_tab/")

gene_order <- as.matrix(row.names(m6A_84sam)) ###
gene_order <- cbind(gene_order,ens[match(gene_order[,1],ens[,1]),2]);gene_order <- as.data.frame(gene_order)

Colon_cancer_1 <- read.table("Colon_cancer_1.txt",header = T, sep = "\t")
alldata <- cbind(gene_order, Colon_cancer_1[match(gene_order[,2],Colon_cancer_1[,1]),9]);rm(Colon_cancer_1)
alldata <- as.data.frame(alldata);b <- as.numeric(alldata[,3] );alldata <- cbind(alldata,b);alldata <- alldata[,-3];rm(b)
colnames(alldata)[1] <-"Gene";colnames(alldata)[2] <-"ens";colnames(alldata)[3] <-"Colon_cancer_1"

for (i in 2:length(list.files())) {
  gene_exp_persam <- read.table(list.files()[i],header = T, sep = "\t")
  alldata <- cbind(alldata, gene_exp_persam[match(alldata[,2],gene_exp_persam[,1]),9])
  colnames( alldata )[i+2]<- as.character( str_split_fixed(c(list.files()[i]), ".txt", 2)[1] ) 
}   ###
rm(gene_exp_persam)

alldata <- na.omit(alldata)###
row.names(alldata)<- alldata[,1]; gene_84sam <- alldata[,-c(1,2)];colnames(gene_84sam)[10] <-"Colon_cancer_10";colnames(gene_84sam)[20] <-"Colon_normal_10"

gene_o <- row.names(gene_84sam)

m6A_84sam<- m6A_84sam[c(gene_o),]
gene_order <- as.data.frame(row.names(m6A_84sam));colnames(gene_order)[1] <- "m6a"

###############
library("preprocessCore")
gene_84sam<- normalize.quantiles(as.matrix(gene_84sam))
colnames(gene_84sam) <- colnames(m6A_84sam);row.names(gene_84sam) <- row.names(m6A_84sam)


p_value<-cor.test(as.numeric(m6A_84sam[1,]),as.numeric(gene_84sam[1,]))[["p.value"]]
rm(i)
for (i in 2:nrow(m6A_84sam)) {
  p <- cor.test(as.numeric(m6A_84sam[i,]),as.numeric(gene_84sam[i,]))[["p.value"]]
  p_value <- c(p_value,p)
}   ###
p_value_adj<- p.adjust(p_value,method ="fdr")

r_value<-cor.test(as.numeric(m6A_84sam[1,]),as.numeric(gene_84sam[1,]))[["estimate"]][["cor"]]
rm(i)
for (i in 2:nrow(m6A_84sam)) {
  r <- cor.test(as.numeric(m6A_84sam[i,]),as.numeric(gene_84sam[i,]))[["estimate"]][["cor"]]
  r_value <- c(r_value,r)
}   ###

fea_m6A_P_R <- cbind(gene_order,p_value,p_value_adj,r_value);colnames(fea_m6A_P_R)[1] <-"Gene"

m6A_loca<- read.csv("./support_files/TCGA_31type_meanTPM5_gene.csv",header = T,row.names = 1)

P_R_gene <- cbind(m6A_loca,fea_m6A_P_R[match(m6A_loca[,1],fea_m6A_P_R[,1]),2]);colnames(P_R_gene)[2] <- "Pvalue"
P_R_gene <- cbind(P_R_gene,fea_m6A_P_R[match(P_R_gene[,1],fea_m6A_P_R[,1]),3]);colnames(P_R_gene)[3] <- "Pvalue_adj"
P_R_gene <- cbind(P_R_gene,fea_m6A_P_R[match(P_R_gene[,1],fea_m6A_P_R[,1]),4]);colnames(P_R_gene)[4] <- "Rvalue"

last_data <- P_R_gene[order(-P_R_gene[,4]),]

library(ggpubr)
ggline(last_data, x = "aaa", y = "Rvalue")


