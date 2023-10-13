library(pheatmap)

heapmap_order <- read.csv("./support_files/Heatmap_data_orderbysamplenum.csv",header = T);heapmap_order<-heapmap_order[,-1]  ###

Sample_Type <- read.csv("./support_files/SampleType.csv",header = T)  ###
annotation_col = data.frame(Type = Sample_Type$Type)#
rownames(annotation_col) = Sample_Type$Sample



pdf(file = "./output_files/heatmap_orderbysamplenum.pdf")

pheatmap(heapmap_order,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = FALSE,cluster_cols = FALSE,
         show_rownames = FALSE,color=colorRampPalette(c("darkgreen","white","darkorange"))(1000))   #

dev.off()


