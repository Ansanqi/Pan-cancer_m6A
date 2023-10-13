library(ggpubr)

m6a_level <- read.csv("./support_files/feature_m6A_level_new.csv",header = T,row.names = 1); m6a_level <- t(m6a_level)
m6a_level2<-data.frame(m6a_level,stringsAsFactors = F);m6a_level2=as.data.frame(lapply(m6a_level2,as.numeric)) #
rownames(m6a_level2) <- rownames(m6a_level)

Sample_Type <- read.csv("./support_files/SampleType3.csv",header = T)  

m6a_boxplot <- cbind(m6a_level2, Sample_Type[match(rownames(m6a_level2),Sample_Type[,1]),c(2,3)])



STAT3 <- m6a_boxplot[,c(37,38,1)] 

p <- ggboxplot(STAT3, x = "Type", y = "STAT3",
               color = "Group", palette = "nejm", 
               add = "none",xlab = FALSE,x.text.angle=45,size =0.9) # 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 5,method = "wilcox.test")    

