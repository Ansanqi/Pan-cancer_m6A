#####peak percent#####

library(gplots)
library(ggplot2)
peak_percent<- read.csv("./support_files/persent_barplot_new.csv",header = T)

pdf(file = "./output_files/percent_Proportion.pdf",width = 5,height = 5)

ze_barplot1 <- barplot(peak_percent[,3] ,beside=T ,main="Cancer Genic Locations" ,legend.text=T,col=c( "darkorange"),  ylab="Proportion",ylim = c(0,0.6),names.arg =peak_percent$Cancer,las = 1 )
text(x = c(0.7,1.9,3.1,4.3),y = 0.02 + peak_percent[,3],labels = peak_percent$Count,adj = c(0.5, 0))

ze_barplot2 <- barplot(peak_percent[,6] ,beside=T ,main="Normal Genic Locations" , legend.text=T,col=c( "darkorange"),  ylab="Proportion",ylim = c(0,0.6),names.arg =peak_percent$Normal,las = 1 )
text(x = c(0.7,1.9,3.1,4.3),y = 0.02 + peak_percent[,6],labels = peak_percent$Count.1,adj = c(0.5, 0))

dev.off()

