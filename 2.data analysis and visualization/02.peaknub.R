#####peak number calculation#####

library(gplots)
library(ggplot2)
library(RColorBrewer)
m6Adata<- read.csv("./support_files/peaknumber.csv",header = T)


pdf(file = "./output_files/peaknumber.pdf",width = 6,height = 6)

ze_barplot <- barplot(m6Adata[,2] ,beside=T , legend.text=T,col=c( "#C59535","#51745E"),  ylab="Peak number",ylim = c(0,16000),names.arg =m6Adata$Group.1,las = 2 )
arrows(ze_barplot,m6Adata$x-m6Adata$x.1*2, ze_barplot, m6Adata$x+m6Adata$x.1*2,length=0.1, angle=90, code=3)

dev.off()
