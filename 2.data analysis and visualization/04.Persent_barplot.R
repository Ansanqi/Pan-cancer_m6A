library(readxl) # ╪сть╟Э
library(ggplot2)
library(tidyr)
library(RColorBrewer)

persent_bar <- read.csv("./support_files/persent.csv",header = T)  
persent_bar2 <- gather(persent_bar, Type, Persent, -Sample)

display.brewer.all()
color <- brewer.pal(12,"Paired")
color <- c("#A6CEE3","#FDBF6F","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#1F78B4","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")

pdf(file = "./output_files/persent_barplot.pdf",width = 10,height = 5)
ggplot(persent_bar2) +
  geom_bar(aes(x = Sample, y = Persent, fill = Type,width=0.6),
           stat = "identity",position = 'fill')+
  scale_fill_manual(values = color) +
  labs(x = "Sample", y = NULL, fill = "Gene Type")+
  coord_flip()
dev.off()

