library('stringr')

m6A_loca<- read.csv("./support_files/m6A_sites_84sample.csv",header = T)

gene_loca <- read.csv("./support_files/gene_RBP_new_new.csv",header = T)

Sample_order <- read.table("./support_files/Sample_order.txt",header = F)
Sample_order <- Sample_order[,1]

############
m6A_forRBP <- cbind(m6A_loca[,1],m6A_loca[,Sample_order])

###########
gene <- as.data.frame(gene_loca[,c(1,3)]) ;colnames(gene)[1] <- "gene";colnames(gene)[2] <- "ens"
setwd("./support_files/gene_abund_tab/")
Colon_cancer_1 <- read.table("Colon_cancer_1.txt",header = T, sep = "\t")
gene_forRBP <- cbind(gene, Colon_cancer_1[match(gene[,2],Colon_cancer_1[,1]),9]);rm(Colon_cancer_1)
colnames(gene_forRBP)[1] <-"Gene";colnames(gene_forRBP)[2] <-"ens";colnames(gene_forRBP)[3] <-"Colon_cancer_1"

for (i in 2:length(list.files())) {
  gene_exp_persam <- read.table(list.files()[i],header = T, sep = "\t")
  gene_forRBP <- cbind(gene_forRBP, gene_exp_persam[match(gene_forRBP[,2],gene_exp_persam[,1]),9])
  colnames( gene_forRBP )[i+2]<- as.character( str_split_fixed(c(list.files()[i]), ".txt", 2)[1] ) 
}   ###
rm(gene_exp_persam)
row.names(gene_forRBP)<- gene_forRBP[,1]; gene_forRBP <- gene_forRBP[,-c(1,2)];colnames(gene_forRBP)[10] <-"Colon_cancer_10";colnames(gene_forRBP)[20] <-"Colon_normal_10"
gene_forRBP <- gene_forRBP[,Sample_order]

write.csv(m6A_forRBP,file = "./output_files/m6a_forRBP.csv")
write.csv(gene_forRBP,file = "./output_files/gene_forRBP.csv")





