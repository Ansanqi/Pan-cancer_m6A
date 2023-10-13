corrdata<- read.csv("./support_files/correlation.csv",header = T,row.names = 1)

corr<- cor(corrdata,use="pairwise.complete.obs",method = "pearson")

library(corrplot)
library(RColorBrewer)
display.brewer.all()

col <- colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))

pdf(file = "./output_files/corralation_new.pdf",width = 10,height = 10)
corrplot(corr, method = "shade",order = "original",tl.cex=0.6,tl.col='black',col = col(100))
dev.off()
