library('stringr')

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

##############

p_value<-cor.test(as.numeric(m6A_84sam[1,]),as.numeric(gene_84sam[1,]))[["p.value"]]
rm(i)
for (i in 2:nrow(m6A_84sam)) {
  p <- cor.test(as.numeric(m6A_84sam[i,]),as.numeric(gene_84sam[i,]))[["p.value"]]
  p_value <- c(p_value,p)
}   ###

p_value_adj<- p.adjust(p_value,method ="fdr")

fea_m6A_P <- cbind(gene_order,p_value,p_value_adj);colnames(fea_m6A_P)[1] <-"Gene"

adj_gene0.01<-subset(fea_m6A_P,fea_m6A_P[,3] <0.01)

write.csv(adj_gene0.01,file = "./output_files/afterFDRadj_m6A2649.csv")
write.csv(fea_m6A_P,file = "./output_files/allm6a_pvalue.csv")


gene<-adj_gene0.01[,1]

Sample_o <- read.table("./support_files/Sample_order.txt",header = F);Sample_o <- Sample_o[,1]
gene_2649 <- gene_84sam[c(gene),c(Sample_o)]
m6A_2649 <- m6A_84sam[c(gene),c(Sample_o)]

gene_2649 <- na.omit(gene_2649)
m6A_2649 <- na.omit(m6A_2649)

write.csv(m6A_2649,file = "./output_files/m6A2649forheatmap.csv")
write.csv(gene_2649,file = "./output_files/gene2649forheatmap.csv")






