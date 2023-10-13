
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(dplyr) 
library(RColorBrewer)
exp <- read.csv("./support_files/RBP_new_TCGA_31type_exp.csv",header = T,row.names = 1)
exp <- na.omit(exp)

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

exp_cv <- as.data.frame(apply(exp, 1, cal_cv)) ;colnames(exp_cv) <- "CV"
exp_cv <- as.matrix(exp_cv)

pdf(file = "./output_files/CV barplot .pdf",height = 3,width = 6)
barplot(exp_cv,horiz = F,beside=T ,ylim = c(0,2), ylab="CV",space = rep(0.2,nrow(exp)))
dev.off()
#
#boxplot
cv_box <- as.data.frame(exp_cv)
gene_type <-  read.csv("./support_files/gene_RBP.csv",header = T)

cv_box <- cbind(cv_box,gene_type[match(row.names(cv_box),gene_type[,1]),2])
colnames(cv_box) <- c("CV_Value","Gene_Type")



pdf(file = "./output_files/boxplot_CV.pdf",height = 5,width = 3.2)
boxplot( CV_Value ~ Gene_Type, cv_box,col = c("#3CB371", "#FF6347"))
dev.off()

a <- cv_box[c(1:6),1];b <- cv_box[c(7:21),1]

wilcox.test(a,b)

library( ggstatsplot)

pdf(file = "./output_files/boxplot_CV_new.pdf",height = 5,width = 5)
ggbetweenstats( cv_box,
                Gene_Type, CV_Value, outlier.tagging=F)
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
# Ìí¼Ó»Ø¹éÏß
abline(lm(y ~ x), col = "black",lwd=1.5)
legend('topleft', cex = 0.5,
       legend = c(paste("P-value =",P),paste("R-value =" ,R)))
dev.off()




