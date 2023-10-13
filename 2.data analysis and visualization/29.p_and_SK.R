

###
m6A <- read.csv("./support_files/m6a_forRBP.csv",header = T,row.names = 1)
gene <- read.csv("./support_files/gene_forRBP.csv",header = T,row.names = 1)

#R_pearson <- cor(m6A,gene,method = 'pearson')

library(Hmisc) 
rcorr_test_pearson <- rcorr(as.matrix(gene),as.matrix(m6A),type = 'pearson')#
rcorr_test_pearson_R <- rcorr_test_pearson$r[1:37,38:2133]
rcorr_test_pearson_P <- rcorr_test_pearson$P[1:37,38:2133]#
write.csv(rcorr_test_pearson_P,'./output_files/P.csv')


##########################

sangji <- read.csv("./output_files/SK.csv",header = T)

library(reshape2)
sangji_new <- melt(sangji, id.vars = "X")

colnames(sangji_new) <- c("target","source","value")

write.csv(sangji_new,file = "./SK2.csv")




