library("pheatmap")
P <- read.csv("./output_files/p.csv",header = T,row.names = 1)
R <- read.csv("./cor_RBP.csv",header = T,row.names = 1)

P <- P[row.names(R),]

#


R[P > 0.05] <- 0  



gene_type <- read.csv("gene_RBP.csv",header = T,row.names = 1);  ###

annotation_row = data.frame(Description = gene_type$Description)#
rownames(annotation_row) = rownames(gene_type)

ann_colors<- list(Description = c(classical = "#3CB371", referenced = "#FF6347"))

bk <- c(seq(-1,-0.0001,by=0.0001),seq(0,1,by=0.0001))

pdf(file = "./cor_RBP_showgenename_new.pdf")
pheatmap(R,border = F,show_colnames = FALSE,annotation_row = annotation_row,cluster_rows = F,cluster_cols = T,annotation_colors = ann_colors,
         show_rownames = T,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-1,1,0.2),
         breaks=bk)   #
dev.off()

