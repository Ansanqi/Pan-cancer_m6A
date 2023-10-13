library(sqldf)
library(tidyr)
library(pheatmap)

###############K=6############
component <- read.csv("./support_files/survival_mat_6con.csv",header = T,row.names = 1)
zdata=sqldf("select Cnew,Type,count(1) from component group by Cnew,Type")
zdata <- na.omit(zdata)
colnames(zdata)[3] <- "count"
dat1=spread(zdata,key=Type,value=count)
dat1[is.na(dat1)] <- 0
row.names(dat1) <- dat1[,1];dat1 <- dat1[,-1]
dat2=matrix(rep(colSums(dat1),nrow(dat1)),nrow = 6,byrow = T)
dat3 <- dat1/dat2
bk = unique(c(seq(0,1, length=1000))) ####
pdf(file = "./output_files/comp_K6.pdf")
pheatmap(dat3,border = T,show_colnames = T,cluster_rows = F,cluster_cols = F,
         show_rownames = T,color=colorRampPalette(c("white","orange","black"))(1000))   #
dev.off()

