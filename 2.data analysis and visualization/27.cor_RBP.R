library(pheatmap)
library(RColorBrewer)

RBP_cor <- read.csv("./support_files/cor_RBP.csv",header = T,row.names = 1);  ###
gene_type <- read.csv("./support_files/gene_RBP_new_new.csv",header = T,row.names = 1);  ###

annotation_row = data.frame(Description = gene_type$Description)#
rownames(annotation_row) = rownames(gene_type)

ann_colors<- list(Description = c(classical = "#3CB371", referenced = "#FF6347"))

pdf(file = "./output_files/cor_RBP_showgenename.pdf")
pheatmap(RBP_cor,border = F,show_colnames = FALSE,annotation_row = annotation_row,cluster_rows = F,cluster_cols = T,annotation_colors = ann_colors,
         show_rownames = T,color=colorRampPalette(c("blue","white","red"))(1000))   #
dev.off()

row_mean = as.data.frame(apply(RBP_cor,1,mean)) ;colnames(row_mean) <- "mean"

row_mean$abs <- abs(row_mean[,1])
row_mean$num <- c(1:nrow(row_mean))
row_mean <- cbind(row_mean,gene_type[match(row.names(row_mean),row.names(gene_type)),1]);colnames(row_mean)[4] <- "type"

pdf(file = "./output_files/abs.pdf")
ggline(row_mean, x = "num", y = "abs",line.size =10 ,point.size=0.5,color = "type",palette = c(classical = "#3CB371", referenced = "#FF6347"))
dev.off()





