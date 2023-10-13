library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(dplyr) 
library(RColorBrewer)
cor <- read.csv("./support_files/cor_RBP.csv",header = T,row.names = 1)
cor <- na.omit(cor)

row_mean = apply(abs(cor),1,mean) ###


pdf(file = "./output_files/barplot.pdf",height = 3,width = 6)
barplot(row_mean,horiz = F,beside=T ,ylim = c(0,0.3), ylab="Absolute value of Correlation")
dev.off()


m6A_express <- log2(m6A_express+1)
gene <- read.csv("./support_files/gene_forRBP.csv",header = T,row.names = 1)
gene <- log2(gene+1)
m6A_express <- t(m6A_express);gene <- t(gene)

row_mean = apply(m6A_express,1,mean)

##############CAPRIN1##############
CAPRIN1 <- cbind(row_mean,gene[,c("CAPRIN1")])

x <- CAPRIN1[,1]
y <- CAPRIN1[,2]
P <- cor.test(x,y)[["p.value"]]
R <- cor.test(x,y)[["estimate"]]

pdf(file = "./output_files/CAPRIN1.pdf",height = 5,width = 5)
plot(x, y, main = "correlation of CAPRIN1 expression and m6A level",
     xlab = "m6A-express level", ylab = "TPM of CAPRIN1",
     pch = 1,cex=0.5, frame = T,col ="#000080")
# 
abline(lm(y ~ x), col = "black",lwd=1.5)
legend('topleft', cex = 0.5,
       legend = c(paste("P-value =",P),paste("R-value =" ,R)))
dev.off()




