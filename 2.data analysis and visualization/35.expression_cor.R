
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(dplyr) 
library(RColorBrewer)
m6A_level <- read.csv("./support_files/m6a_forRBP.csv",header = T)
m6A_level <- t(m6A_level)
colnames( m6A_level) <- m6A_level[1,];m6A_level <- m6A_level[-1,]
m6A_level <- as.data.frame(m6A_level);

for (i in 1:ncol(m6A_level)) {
  m6A_level[,i] <- as.numeric(m6A_level[,i])
}  

row_mean = apply(m6A_level,1,mean)

gene <- read.csv("./support_files/gene_forRBP.csv",header = T,row.names = 1)
gene <- log2(gene+1)
gene <- t(gene)


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

