


library(ggplot2)
#KEGG 
pathway = read.csv("./support_files/KEGG.csv",header = T)
library(ggplot2) 

p = ggplot(pathway,aes(-1*log10(PValue),Term2))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(PValue)))
pr = pbubble+scale_color_gradient(low="green",high = "red")
pr = pr+labs(color=expression(-log[10](PValue)),size="Gene number",x="PValue",y="Pathway name",title="KEGG Pathway enrichment")+theme(plot.title=element_text(hjust=0.5))
pr=pr + theme_bw()
pdf(file = "./output_files/KEGG.pdf",height = 6,width = 7)
pr
dev.off()



#GO
pathway = read.csv("./support_files/BP.csv",header = T)
library(ggplot2) 

p = ggplot(pathway,aes(-1*log10(PValue),Term2))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(PValue)))
pr = pbubble+scale_color_gradient(low="green",high = "red")
pr = pr+labs(color=expression(-log[10](PValue)),size="Gene number",x="PValue",y="BP Term",title="GO BP enrichment")+theme(plot.title=element_text(hjust=0.5))
pr=pr + theme_bw()
pdf(file = "./output_files/BP.pdf",height = 6,width = 8.4)
pr
dev.off()
################################################

pathway = read.csv("./support_files/CC.csv",header = T)
library(ggplot2) 

p = ggplot(pathway,aes(-1*log10(PValue),Term2))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(PValue)))
pr = pbubble+scale_color_gradient(low="green",high = "red")
pr = pr+labs(color=expression(-log[10](PValue)),size="Gene number",x="PValue",y="CC Term",title="GO CC enrichment")+theme(plot.title=element_text(hjust=0.5))
pr=pr + theme_bw()
pdf(file = "./output_files/CC.pdf",height = 6,width = 7)
pr
dev.off()


###############################################

pathway = read.csv("./support_files/MF.csv",header = T)
library(ggplot2) 

p = ggplot(pathway,aes(-1*log10(PValue),Term2))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(PValue)))
pr = pbubble+scale_color_gradient(low="green",high = "red")
pr = pr+labs(color=expression(-log[10](PValue)),size="Gene number",x="PValue",y="MF Term",title="GO MF enrichment")+theme(plot.title=element_text(hjust=0.5))
pr=pr + theme_bw()
pdf(file = "./output_files/MF.pdf",height = 6,width = 10.5)
pr
dev.off()


