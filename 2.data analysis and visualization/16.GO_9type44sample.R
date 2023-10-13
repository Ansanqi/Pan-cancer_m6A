
library(ggplot2)
library(RColorBrewer)

GO_BP <- read.csv("./support_files/GO_BP_9type.csv",header = T)
GO_CC <- read.csv("./support_files/GO_CC_9type.csv",header = T)
GO_MF <- read.csv("./support_files/GO_MF_9type.csv",header = T)

display.brewer.all()
color <- c("#F08D85", "#D684B4", "#A8BD1C", "#C894C1", "#85AFDD","#E7A815", "#4EB45B", "#45BFE5","#4DBAA7")
colorl <- rep(color,times = c(5,10,10,6,10,10,10,8,8))
colorl

ggplot(GO_BP) +
  aes(x = name2, y = X.log10Pvalue, fill = Type) +
  geom_bar(stat = "identity",colour="black") +
  #scale_fill_hue() +
  coord_flip()+
  scale_fill_manual(values =color)+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.direction = "vertical",
    legend.position = c(0.9,0.92),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "black"),
    plot.background = element_blank(),
    
  )#+geom_hline(aes(yintercept=-log10(0.05)),linetype=5,col="red")
dev.off()

