library(GSVA)
library(GSEABase)
library(pheatmap)


c2gmt <- getGmt("./support_files/c2.cp.kegg.v7.5.1.symbols.gmt")

m6A_Pan <- read.csv("./support_files/ssGSEA.csv" , header = T, row.names = 1) ###

gs.exp <- gsva(as.matrix(m6A_Pan), c2gmt, kcdf = "Poisson", min.sz = 3)

write.csv(gs.exp, file = "./9cancer_KEGG_score.csv")

gs.exp_order <- read.csv("./9cancer_KEGG_score_order.csv" , header = T, row.names = 1) 
gs.exp_order <- gs.exp_order[,-c(45:47)]

Sample_Type <- read.csv("./support_files/SampleType.csv",header = T,row.names = 1)  ###

bk = unique(c(seq(-1.5,1.5, length=100))) ####

pheatmap(gs.exp_order,
         show_colnames = F,scale = "row",border = F,
         cluster_rows = F,cluster_cols = F,
         breaks=bk,annotation_col = Sample_Type,
         treeheight_row = 20,cellheight=6,cellwidth=6,
         fontsize=6,gaps_col = c(10,14,17,19,23,31,37,42),width = 12)
dev.off()


