library('stringr')

m6A_forRBP<- read.csv("./output_files/m6a_forRBP.csv",header = T)
gene_forRBP<- read.csv("./output_files/gene_forRBP.csv",header = T,row.names = 1)

sam <- colnames(gene_forRBP)

random1 <- sample(sam,length(sam));random1 <- c("X",random1)
random2 <- sample(sam,length(sam))

m6A_random <- m6A_forRBP[,random1]
gene_random <- gene_forRBP[,random2]


write.csv(m6A_random,file = "./output_files/m6A_randomsample_1.csv")
write.csv(gene_random,file = "./output_files/gene_randomsample_1.csv")


random1 <- sample(sam,length(sam));random1 <- c("X",random1)
random2 <- sample(sam,length(sam))

m6A_random <- m6A_forRBP[,random1]
gene_random <- gene_forRBP[,random2]

write.csv(m6A_random,file = "./output_files/m6A_randomsample_2.csv")
write.csv(gene_random,file = "./output_files/gene_randomsample_2.csv")

random1 <- sample(sam,length(sam));random1 <- c("X",random1)
random2 <- sample(sam,length(sam))

m6A_random <- m6A_forRBP[,random1]
gene_random <- gene_forRBP[,random2]

write.csv(m6A_random,file = "./output_files/m6A_randomsample_3.csv")
write.csv(gene_random,file = "./output_files/gene_randomsample_3.csv")
