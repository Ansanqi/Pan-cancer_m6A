library('stringr')

m6A_loca<- read.csv("./support_files/m6A_sites.csv",header = T)

gene_loca <- read.csv("./support_files/gene_RBP.csv",header = T)

Sample_order <- read.table("./support_files/Sample_order.txt",header = F)
Sample_order <- Sample_order[,1]

m6A_forRBP <- cbind(m6A_loca[,1],m6A_loca[,Sample_order])
m6A_forRBP <- t(m6A_forRBP)
colnames(m6A_forRBP) <- m6A_forRBP[1,];m6A_forRBP <- t(m6A_forRBP)
m6A_forRBP <- m6A_forRBP[,-1];m6A_forRBP <- as.data.frame(m6A_forRBP)
m6A_forRBP2 <- as.data.frame(lapply(m6A_forRBP,as.numeric));row.names(m6A_forRBP2) <- row.names(m6A_forRBP);rm(m6A_forRBP)
m6A_forRBP <- m6A_forRBP2;rm(m6A_forRBP2)


gene <- as.data.frame(gene_loca[,c(1,3)]) ;colnames(gene)[1] <- "gene";colnames(gene)[2] <- "ens"
setwd("./gene_abund_tab/")
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

##
m6A_forRBP <- t(m6A_forRBP[,Sample_order]);gene_forRBP <- t(gene_forRBP[,Sample_order])
library(Hmisc) 
rcorr_test_pearson <- rcorr(as.matrix(gene_forRBP),as.matrix(m6A_forRBP),type = 'pearson')
R <- rcorr_test_pearson$r[1:52,53:2148]
P <- rcorr_test_pearson$P[1:52,53:2148]
rm(rcorr_test_pearson)
##
a <- R[c(1:20),];row_mean = apply(a,1,mean);a <- cbind(a,row_mean);a <- a[order(a[,ncol(a)],decreasing=T),];a <- a[,-ncol(a)]
b <- R[c(21:nrow(R)),];row_mean = apply(b,1,mean);b <- cbind(b,row_mean);b <- b[order(b[,ncol(b)],decreasing=T),];b <- b[,-ncol(b)];rm(row_mean)
R <- rbind(a,b);rm(a,b)

P <- P[row.names(R),colnames(R)]
R[P > 0.05] <- 0  

gene_type <- read.csv("./support_files/gene_RBP.csv",header = T,row.names = 1);  ###

annotation_row = data.frame(Description = gene_type$Description)#
rownames(annotation_row) = rownames(gene_type)

ann_colors<- list(Description = c(classical = "#3CB371", referenced = "#FF6347"))

bk <- c(seq(-1,-0.0001,by=0.0001),seq(0,1,by=0.0001))
library(pheatmap)
pdf(file = "./output_files/cor_RBP_showgenename_new_32referenced.pdf")
pheatmap(R,border = F,show_colnames = FALSE,annotation_row = annotation_row,cluster_rows = F,cluster_cols = T,annotation_colors = ann_colors,
         show_rownames = T,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-1,1,0.2),
         breaks=bk)   #
dev.off()

###########CDF
R2 <- cor(gene_forRBP,m6A_forRBP,method = 'pearson')

#
random1 <- sample(Sample_order,length(Sample_order))
random2 <- sample(Sample_order,length(Sample_order))

m6A_random <- m6A_forRBP[random1,]
gene_random <- gene_forRBP[random2,]
R3 <- cor(gene_random,m6A_random,method = 'pearson')

pdf(file = "./output_files/CDF.pdf")
plot( ecdf( abs(as.numeric( as.matrix( R2[c(1:20),] )) ))    , col="orange", xlim=c(0.01,0.9)
      ,cex=0.001,xlab="Absolute value of the correlation between RBP expression and cancer specific m6A level" , ylab = "Cumulative Density Function")
plot( ecdf( abs(as.numeric( as.matrix(R2[c(21:52),] ))))  ,add=T,col="darkred"   )
plot( ecdf( abs(as.numeric( as.matrix(R3 ))))  ,add=T,col="black"   )
dev.off()

########
R4 <- R2[row.names(R),]
R4 <- na.omit(R4)
row_mean = apply(abs(R4),1,mean) ###

pdf(file = "./output_files/barplot.pdf",height = 3,width = 6)
barplot(row_mean,horiz = F,beside=T ,ylim = c(0,0.3), ylab="Absolute value of Correlation Coefficient")
dev.off()

#########
P2 <- t(P)
P2[P2 > 0.05] <- NA  
write.csv(P2,'./output_files/P.csv')


sangji <- read.csv("./output_files/SK.csv",header = T)

library(reshape2)
sangji_new <- melt(sangji, id.vars = "X")

colnames(sangji_new) <- c("target","source","value")

write.csv(sangji_new,file = "./SK2.csv")

