library('stringr')
library("pheatmap")
library("dplyr")

m6A_84_alldata<- read.csv("./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)
sample_order <- read.table("./support_files/SampleType.csv",header = T,sep = ",")
m6A_84 <- m6A_84_alldata[,-c(2:8)]
m6A_84[is.na(m6A_84)] <- 0
m6A_84 <- unique(m6A_84)#
m6A_84 <- m6A_84 %>% distinct(rownames.output., .keep_all = T)

row.names(m6A_84) <- m6A_84[,1]
m6A_84 <- m6A_84[,-1]
m6A_84 <- m6A_84[,sample_order[,1]]

cancer_m6A <- m6A_84[,c(1:45)]
normal_m6A <- m6A_84[,c(46:84)]


p_value<-t.test(as.numeric(cancer_m6A[1,]),as.numeric(normal_m6A[1,]))[["p.value"]]

for (i in 2:nrow(m6A_84)) {
  p <- t.test(as.numeric(cancer_m6A[i,]),as.numeric(normal_m6A[i,]))[["p.value"]]
  p_value <- c(p_value,p)
}   ###

p_value_adj<- p.adjust(p_value,method ="fdr")
fea_m6A_P <- cbind(p_value,p_value_adj);row.names(fea_m6A_P) <- row.names(m6A_84)
adj_sites0.05<-subset(fea_m6A_P,fea_m6A_P[,2] <0.05)

sites <- row.names(adj_sites0.05)#

m6A_84_4map <- m6A_84[sites,]

ann_colors<- list(Type = c(cancer="#DA4C35",normal="#4BB1C9"))

annotation_col = data.frame(Type = sample_order$Type)#
rownames(annotation_col) = sample_order$Sample

p <- pheatmap(m6A_84_4map,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = TRUE,cluster_cols = F,
              show_rownames = FALSE,color=colorRampPalette(c("darkgreen","white","darkorange"))(1000),annotation_colors = ann_colors)   #

pdf(file = "./output_files/m6A_adj0.05.pdf")
p
dev.off()
dev.off()

sites_in_data <- cbind(c(1:length(sites)),p[["tree_row"]][["order"]]);row.names(sites_in_data) <- sites
sites_in_order <- sites_in_data[order(sites_in_data[,2],decreasing=F),]

pick_gene <- cbind(sites_in_order,m6A_84_alldata[match(row.names(sites_in_order),m6A_84_alldata[,1]),c(3,2)])

gene_order <- unique(pick_gene[,3])

setwd("./gene_abund_tab/")
gene_order <- cbind(gene_order,pick_gene[match(gene_order,pick_gene[,3]),4])

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
row.names(alldata)<- alldata[,1]; gene_84 <- alldata[,-c(1,2)];colnames(gene_84)[10] <-"Colon_cancer_10";colnames(gene_84)[20] <-"Colon_normal_10"


gene_84 <- gene_84[,sample_order[,1]]
gene_84 <- log2(gene_84+1)
library("preprocessCore")
gene_84_4map<- normalize.quantiles(as.matrix(gene_84))
colnames(gene_84_4map) <- colnames(gene_84);row.names(gene_84_4map) <- row.names(gene_84)

pdf(file = "./output_files/gene.pdf")
pheatmap(gene_84_4map,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = F,cluster_cols = F,
         show_rownames = FALSE,color=colorRampPalette(c("darkgreen","white","darkorange"))(1000))   #
dev.off()

write <- cbind(pick_gene,adj_sites0.05[match(row.names(pick_gene),row.names(adj_sites0.05)),c(1,2)])

write.csv(write,file = "./genelist_0.05.csv")



cancer_m <- m6A_84_4map[,c(1:45)]
normal_m <- m6A_84_4map[,c(46:84)]

cal_cv=function(x){  # 
  y=na.omit(x)
  return(sd(y)/mean(y))
}

cancer_cv <- apply(cancer_m, 1, cal_cv)
normal_cv <- apply(normal_m, 1, cal_cv)

wilcox.test(cancer_cv,normal_cv)


pdf(file = "./output_files/CV.pdf",height = 5,width = 3)
boxplot(cancer_cv,normal_cv,col = c("#DA4C35", "#4BB1C9"))
dev.off()



