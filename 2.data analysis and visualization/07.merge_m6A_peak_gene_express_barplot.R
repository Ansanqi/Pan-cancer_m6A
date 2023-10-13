library('stringr')
library(ggplot2)
library(ggpubr)


load("./support_files/peak_mean_max_forplot.RData")

pdf(file = "./output_files/TPM_t.test.pdf",height = 5,width = 3)
p_value1<- t.test(cancer_mean,normal_mean)[["p.value"]]
boxplot(cancer_mean,normal_mean,names = c("cancer","normal"),col = c("#FF4500", "#008080"),boxwex = 0.5,notch = T,ylab="log2 TPM",main=c("t.test=",p_value1))
dev.off()


pdf(file = "./output_files/Max.pdf",height = 5,width = 3)
p_value2<- wilcox.test(cancer_max,normal_max)[["p.value"]]
boxplot(cancer_max,normal_max,names = c("cancer","normal"),col = c("#FF4500", "#008080"),boxwex = 0.5,notch = T,ylab="Maximum m6A ratio",main=c("wilcox.test=",p_value2))
dev.off()
