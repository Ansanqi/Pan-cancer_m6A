gene_order <- read.table("./support_files/Gene_order.txt",header = T) 

Colon_cancer_1 <- read.table("./support_files/gene_abund_tab/SRR14930904gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(gene_order, Colon_cancer_1[match(gene_order[,2],Colon_cancer_1[,1]),9]);rm(Colon_cancer_1);rm(gene_order)
colnames(alldata)[3] <-"Colon_cancer_1"

Colon_cancer_2 <- read.table("./support_files/gene_abund_tab/SRR14930905gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_2[match(alldata[,2],Colon_cancer_2[,1]),9]);rm(Colon_cancer_2)
colnames(alldata)[4] <-"Colon_cancer_2"

Colon_cancer_3 <- read.table("./support_files/gene_abund_tab/SRR14930906gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_3[match(alldata[,2],Colon_cancer_3[,1]),9]);rm(Colon_cancer_3)
colnames(alldata)[5] <-"Colon_cancer_3"

Colon_cancer_4 <- read.table("./support_files/gene_abund_tab/SRR14930907gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_4[match(alldata[,2],Colon_cancer_4[,1]),9]);rm(Colon_cancer_4)
colnames(alldata)[6] <-"Colon_cancer_4"

Colon_cancer_5 <- read.table("./support_files/gene_abund_tab/SRR14930908gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_5[match(alldata[,2],Colon_cancer_5[,1]),9]);rm(Colon_cancer_5)
colnames(alldata)[7] <-"Colon_cancer_5"

Colon_cancer_6 <- read.table("./support_files/gene_abund_tab/SRR14930909gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_6[match(alldata[,2],Colon_cancer_6[,1]),9]);rm(Colon_cancer_6)
colnames(alldata)[8] <-"Colon_cancer_6"

Colon_cancer_7 <- read.table("./support_files/gene_abund_tab/SRR14930910gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_7[match(alldata[,2],Colon_cancer_7[,1]),9]);rm(Colon_cancer_7)
colnames(alldata)[9] <-"Colon_cancer_7"

Colon_cancer_8 <- read.table("./support_files/gene_abund_tab/SRR14930911gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_8[match(alldata[,2],Colon_cancer_8[,1]),9]);rm(Colon_cancer_8)
colnames(alldata)[10] <-"Colon_cancer_8"

Colon_cancer_9 <- read.table("./support_files/gene_abund_tab/SRR14930912gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_9[match(alldata[,2],Colon_cancer_9[,1]),9]);rm(Colon_cancer_9)
colnames(alldata)[11] <-"Colon_cancer_9"

Colon_cancer_10 <- read.table("./support_files/gene_abund_tab/SRR14930913gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_cancer_10[match(alldata[,2],Colon_cancer_10[,1]),9]);rm(Colon_cancer_10)
colnames(alldata)[12] <-"Colon_cancer_10"

Endometrial_cancer_1 <- read.table("./support_files/gene_abund_tab/SRR5194781gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_cancer_1[match(alldata[,2],Endometrial_cancer_1[,1]),9]);rm(Endometrial_cancer_1)
colnames(alldata)[13] <-"Endometrial_cancer_1"

Endometrial_cancer_2 <- read.table("./support_files/gene_abund_tab/SRR5194785gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_cancer_2[match(alldata[,2],Endometrial_cancer_2[,1]),9]);rm(Endometrial_cancer_2)
colnames(alldata)[14] <-"Endometrial_cancer_2"

Endometrial_cancer_3 <- read.table("./support_files/gene_abund_tab/SRR5194789gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_cancer_3[match(alldata[,2],Endometrial_cancer_3[,1]),9]);rm(Endometrial_cancer_3)
colnames(alldata)[15] <-"Endometrial_cancer_3"

Endometrial_cancer_5 <- read.table("./support_files/gene_abund_tab/SRR5194795gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_cancer_5[match(alldata[,2],Endometrial_cancer_5[,1]),9]);rm(Endometrial_cancer_5)
colnames(alldata)[16] <-"Endometrial_cancer_5"

Glioma_1 <- read.table("./support_files/gene_abund_tab/SRR12518107gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Glioma_1[match(alldata[,2],Glioma_1[,1]),9]);rm(Glioma_1)
colnames(alldata)[17] <-"Glioma_1"

Glioma_2 <- read.table("./support_files/gene_abund_tab/SRR12518109gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Glioma_2[match(alldata[,2],Glioma_2[,1]),9]);rm(Glioma_2)
colnames(alldata)[18] <-"Glioma_2"

Glioma_3 <- read.table("./support_files/gene_abund_tab/SRR12518111gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Glioma_3[match(alldata[,2],Glioma_3[,1]),9]);rm(Glioma_3)
colnames(alldata)[19] <-"Glioma_3"

Leukemia_1 <- read.table("./support_files/gene_abund_tab/SRR3066066gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Leukemia_1[match(alldata[,2],Leukemia_1[,1]),9]);rm(Leukemia_1)
colnames(alldata)[20] <-"Leukemia_1"

Leukemia_2 <- read.table("./support_files/gene_abund_tab/SRR3066068gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Leukemia_2[match(alldata[,2],Leukemia_2[,1]),9]);rm(Leukemia_2)
colnames(alldata)[21] <-"Leukemia_2"

Liver_cancer_1 <- read.table("./support_files/gene_abund_tab/SRR7965968gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_cancer_1[match(alldata[,2],Liver_cancer_1[,1]),9]);rm(Liver_cancer_1)
colnames(alldata)[22] <-"Liver_cancer_1"

Liver_cancer_2 <- read.table("./support_files/gene_abund_tab/SRR7965976gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_cancer_2[match(alldata[,2],Liver_cancer_2[,1]),9]);rm(Liver_cancer_2)
colnames(alldata)[23] <-"Liver_cancer_2"

Liver_cancer_3 <- read.table("./support_files/gene_abund_tab/SRR7965988gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_cancer_3[match(alldata[,2],Liver_cancer_3[,1]),9]);rm(Liver_cancer_3)
colnames(alldata)[24] <-"Liver_cancer_3"

Liver_cancer_4 <- read.table("./support_files/gene_abund_tab/SRR7965992gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_cancer_4[match(alldata[,2],Liver_cancer_4[,1]),9]);rm(Liver_cancer_4)
colnames(alldata)[25] <-"Liver_cancer_4"

Lung_cancer_1 <- read.table("./support_files/gene_abund_tab/SRR7535600gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_1[match(alldata[,2],Lung_cancer_1[,1]),9]);rm(Lung_cancer_1)
colnames(alldata)[26] <-"Lung_cancer_1"

Lung_cancer_2 <- read.table("./support_files/gene_abund_tab/SRR7535600gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_2[match(alldata[,2],Lung_cancer_2[,1]),9]);rm(Lung_cancer_2)
colnames(alldata)[27] <-"Lung_cancer_2"

Lung_cancer_3 <- read.table("./support_files/gene_abund_tab/SRR7535603gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_3[match(alldata[,2],Lung_cancer_3[,1]),9]);rm(Lung_cancer_3)
colnames(alldata)[28] <-"Lung_cancer_3"

Lung_cancer_4 <- read.table("./support_files/gene_abund_tab/SRR7535603gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_4[match(alldata[,2],Lung_cancer_4[,1]),9]);rm(Lung_cancer_4)
colnames(alldata)[29] <-"Lung_cancer_4"

Lung_cancer_5 <- read.table("./support_files/gene_abund_tab/SRR7535606gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_5[match(alldata[,2],Lung_cancer_5[,1]),9]);rm(Lung_cancer_5)
colnames(alldata)[30] <-"Lung_cancer_5"

Lung_cancer_6 <- read.table("./support_files/gene_abund_tab/SRR7535606gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_6[match(alldata[,2],Lung_cancer_6[,1]),9]);rm(Lung_cancer_6)
colnames(alldata)[31] <-"Lung_cancer_6"

Lung_cancer_7 <- read.table("./support_files/gene_abund_tab/SRR7535609gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_7[match(alldata[,2],Lung_cancer_7[,1]),9]);rm(Lung_cancer_7)
colnames(alldata)[32] <-"Lung_cancer_7"

Lung_cancer_8 <- read.table("./support_files/gene_abund_tab/SRR7535609gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_cancer_8[match(alldata[,2],Lung_cancer_8[,1]),9]);rm(Lung_cancer_8)
colnames(alldata)[33] <-"Lung_cancer_8"

Ovarian_cancer_1 <- read.table("./support_files/gene_abund_tab/SRR7763558gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_cancer_1[match(alldata[,2],Ovarian_cancer_1[,1]),9]);rm(Ovarian_cancer_1)
colnames(alldata)[34] <-"Ovarian_cancer_1"

Ovarian_cancer_2 <- read.table("./support_files/gene_abund_tab/SRR7763559gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_cancer_2[match(alldata[,2],Ovarian_cancer_2[,1]),9]);rm(Ovarian_cancer_2)
colnames(alldata)[35] <-"Ovarian_cancer_2"

Ovarian_cancer_3 <- read.table("./support_files/gene_abund_tab/SRR7763560gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_cancer_3[match(alldata[,2],Ovarian_cancer_3[,1]),9]);rm(Ovarian_cancer_3)
colnames(alldata)[36] <-"Ovarian_cancer_3"

Ovarian_cancer_4 <- read.table("./support_files/gene_abund_tab/SRR7763561gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_cancer_4[match(alldata[,2],Ovarian_cancer_4[,1]),9]);rm(Ovarian_cancer_4)
colnames(alldata)[37] <-"Ovarian_cancer_4"

Ovarian_cancer_5 <- read.table("./support_files/gene_abund_tab/SRR7763562gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_cancer_5[match(alldata[,2],Ovarian_cancer_5[,1]),9]);rm(Ovarian_cancer_5)
colnames(alldata)[38] <-"Ovarian_cancer_5"

Ovarian_cancer_6 <- read.table("./support_files/gene_abund_tab/SRR7763563gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_cancer_6[match(alldata[,2],Ovarian_cancer_6[,1]),9]);rm(Ovarian_cancer_6)
colnames(alldata)[39] <-"Ovarian_cancer_6"

Salivary_gland_can_1 <- read.table("./support_files/gene_abund_tab/SRR13091695gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_can_1[match(alldata[,2],Salivary_gland_can_1[,1]),9]);rm(Salivary_gland_can_1)
colnames(alldata)[40] <-"Salivary_gland_can_1"

Salivary_gland_can_2 <- read.table("./support_files/gene_abund_tab/SRR13091699gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_can_2[match(alldata[,2],Salivary_gland_can_2[,1]),9]);rm(Salivary_gland_can_2)
colnames(alldata)[41] <-"Salivary_gland_can_2"

Salivary_gland_can_3 <- read.table("./support_files/gene_abund_tab/SRR13091703gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_can_3[match(alldata[,2],Salivary_gland_can_3[,1]),9]);rm(Salivary_gland_can_3)
colnames(alldata)[42] <-"Salivary_gland_can_3"

Salivary_gland_can_4 <- read.table("./support_files/gene_abund_tab/SRR13091707gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_can_4[match(alldata[,2],Salivary_gland_can_4[,1]),9]);rm(Salivary_gland_can_4)
colnames(alldata)[43] <-"Salivary_gland_can_4"

Salivary_gland_can_5 <- read.table("./support_files/gene_abund_tab/SRR13091711gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_can_5[match(alldata[,2],Salivary_gland_can_5[,1]),9]);rm(Salivary_gland_can_5)
colnames(alldata)[44] <-"Salivary_gland_can_5"

Stomach_cancer_1 <- read.table("./support_files/gene_abund_tab/c1_inputgene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Stomach_cancer_1[match(alldata[,2],Stomach_cancer_1[,1]),9]);rm(Stomach_cancer_1)
colnames(alldata)[45] <-"Stomach_cancer_1"

Stomach_cancer_2 <- read.table("./support_files/gene_abund_tab/c3_inputgene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Stomach_cancer_2[match(alldata[,2],Stomach_cancer_2[,1]),9]);rm(Stomach_cancer_2)
colnames(alldata)[46] <-"Stomach_cancer_2"

alldata <- alldata[,-c(1:2)] 

heapmap_gene  <- log(alldata+1,2) 

write.csv(heapmap_gene,file = "./order.csv")

heapmap_gene_ok<- read.csv("./order_ok.csv",header = T)


library(pheatmap)

Sample_Type <- read.csv("./support_files/SampleType.csv",header = T)  ###
annotation_col = data.frame(Type = Sample_Type$Type)#
rownames(annotation_col) = Sample_Type$Sample

library("preprocessCore")
table_dup<- normalize.quantiles(as.matrix(heapmap_gene_ok[,-c(25,27,29,31)]))
library(gplots)
library("RColorBrewer")
palette =rev( colorRampPalette((brewer.pal(n = 11,	name="RdBu")))(300) )
table_dup<- scale(  heapmap_gene_ok[,-c(25,27,29,31)], center = T)

pdf(file = "./output_files/gene_heatmap_new.pdf")

pheatmap(as.matrix(table_dup) ,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = FALSE,cluster_cols = FALSE,
         show_rownames = FALSE, breaks = c(seq(-4,-1.3,length=100),seq(-1.29,1.4,length=100),seq(1.49,4,length=101)), color=palette)   #

dev.off()


