library(ggpubr)

gene_order <- read.table("./support_files/gene_list.txt",header = T) 

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

Lung_normal_4 <- read.table("./support_files/gene_abund_tab/CRR042297gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_normal_4[match(alldata[,2],Lung_normal_4[,1]),9]);rm(Lung_normal_4)
colnames(alldata)[47] <-"Lung_normal_4"

Lung_normal_3 <- read.table("./support_files/gene_abund_tab/CRR055534gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_normal_3[match(alldata[,2],Lung_normal_3[,1]),9]);rm(Lung_normal_3)
colnames(alldata)[48] <-"Lung_normal_3"

Lung_normal_2<- read.table("./support_files/gene_abund_tab/CRR073019gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_normal_2[match(alldata[,2],Lung_normal_2[,1]),9]);rm(Lung_normal_2)
colnames(alldata)[49] <-"Lung_normal_2"

Lung_normal_1 <- read.table("./support_files/gene_abund_tab/CRR073021gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Lung_normal_1[match(alldata[,2],Lung_normal_1[,1]),9]);rm(Lung_normal_1)
colnames(alldata)[50] <-"Lung_normal_1"

Stomach_normal_1 <- read.table("./support_files/gene_abund_tab/n1_inputgene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Stomach_normal_1[match(alldata[,2],Stomach_normal_1[,1]),9]);rm(Stomach_normal_1)
colnames(alldata)[51] <-"Stomach_normal_1"

Stomach_normal_2 <- read.table("./support_files/gene_abund_tab/n2_inputgene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Stomach_normal_2[match(alldata[,2],Stomach_normal_2[,1]),9]);rm(Stomach_normal_2)
colnames(alldata)[52] <-"Stomach_normal_2"

PBMC_1<- read.table("./support_files/gene_abund_tab/PBMC_6_inputgene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, PBMC_1[match(alldata[,2],PBMC_1[,1]),9]);rm(PBMC_1)
colnames(alldata)[53] <-"PBMC_1"

PBMC_2 <- read.table("./support_files/gene_abund_tab/PBMC_7_inputgene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, PBMC_2[match(alldata[,2],PBMC_2[,1]),9]);rm(PBMC_2)
colnames(alldata)[54] <-"PBMC_2"

Endometrial_normal_1 <- read.table("./support_files/gene_abund_tab/SRR5194779gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_normal_1[match(alldata[,2],Endometrial_normal_1[,1]),9]);rm(Endometrial_normal_1)
colnames(alldata)[55] <-"Endometrial_normal_1"

Endometrial_normal_2 <- read.table("./support_files/gene_abund_tab/SRR5194783gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_normal_2[match(alldata[,2],Endometrial_normal_2[,1]),9]);rm(Endometrial_normal_2)
colnames(alldata)[56] <-"Endometrial_normal_2"

Endometrial_normal_3 <- read.table("./support_files/gene_abund_tab/SRR5194787gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_normal_3[match(alldata[,2],Endometrial_normal_3[,1]),9]);rm(Endometrial_normal_3)
colnames(alldata)[57] <-"Endometrial_normal_3"

Endometrial_cancer_4 <- read.table("./support_files/gene_abund_tab/SRR5194791gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_cancer_4[match(alldata[,2],Endometrial_cancer_4[,1]),9]);rm(Endometrial_cancer_4)
colnames(alldata)[58] <-"Endometrial_cancer_4"

Endometrial_normal_4 <- read.table("./support_files/gene_abund_tab/SRR5194793gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_normal_4[match(alldata[,2],Endometrial_normal_4[,1]),9]);rm(Endometrial_normal_4)
colnames(alldata)[59] <-"Endometrial_normal_4"

Endometrial_normal_5 <- read.table("./support_files/gene_abund_tab/SRR5194797gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Endometrial_normal_5[match(alldata[,2],Endometrial_normal_5[,1]),9]);rm(Endometrial_normal_5)
colnames(alldata)[60] <-"Endometrial_normal_5"

Ovarian_normal_1 <- read.table("./support_files/gene_abund_tab/SRR7763564gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_1[match(alldata[,2],Ovarian_normal_1[,1]),9]);rm(Ovarian_normal_1)
colnames(alldata)[61] <-"Ovarian_normal_1"

Ovarian_normal_2 <- read.table("./support_files/gene_abund_tab/SRR7763565gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_2[match(alldata[,2],Ovarian_normal_2[,1]),9]);rm(Ovarian_normal_2)
colnames(alldata)[62] <-"Ovarian_normal_2"

Ovarian_normal_3 <- read.table("./support_files/gene_abund_tab/SRR7763566gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_3[match(alldata[,2],Ovarian_normal_3[,1]),9]);rm(Ovarian_normal_3)
colnames(alldata)[63] <-"Ovarian_normal_3"

Ovarian_normal_4 <- read.table("./support_files/gene_abund_tab/SRR7763567gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_4[match(alldata[,2],Ovarian_normal_4[,1]),9]);rm(Ovarian_normal_4)
colnames(alldata)[64] <-"Ovarian_normal_4"

Ovarian_normal_5 <- read.table("./support_files/gene_abund_tab/SRR7763568gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_5[match(alldata[,2],Ovarian_normal_5[,1]),9]);rm(Ovarian_normal_5)
colnames(alldata)[65] <-"Ovarian_normal_5"

Ovarian_normal_6 <- read.table("./support_files/gene_abund_tab/SRR7763569gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_6[match(alldata[,2],Ovarian_normal_6[,1]),9]);rm(Ovarian_normal_6)
colnames(alldata)[66] <-"Ovarian_normal_6"

Ovarian_normal_7 <- read.table("./support_files/gene_abund_tab/SRR7763570gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Ovarian_normal_7[match(alldata[,2],Ovarian_normal_7[,1]),9]);rm(Ovarian_normal_7)
colnames(alldata)[67] <-"Ovarian_normal_7"

Liver_normal_1 <- read.table("./support_files/gene_abund_tab/SRR7965970gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_normal_1[match(alldata[,2],Liver_normal_1[,1]),9]);rm(Liver_normal_1)
colnames(alldata)[68] <-"Liver_normal_1"

Liver_normal_2 <- read.table("./support_files/gene_abund_tab/SRR7965974gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_normal_2[match(alldata[,2],Liver_normal_2[,1]),9]);rm(Liver_normal_2)
colnames(alldata)[69] <-"Liver_normal_2"

Liver_normal_3 <- read.table("./support_files/gene_abund_tab/SRR7965978gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_normal_3[match(alldata[,2],Liver_normal_3[,1]),9]);rm(Liver_normal_3)
colnames(alldata)[70] <-"Liver_normal_3"

Liver_normal_4 <- read.table("./support_files/gene_abund_tab/SRR7965994gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Liver_normal_4[match(alldata[,2],Liver_normal_4[,1]),9]);rm(Liver_normal_4)
colnames(alldata)[71] <-"Liver_normal_4"

Glioma_nor_1 <- read.table("./support_files/gene_abund_tab/SRR12518113gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Glioma_nor_1[match(alldata[,2],Glioma_nor_1[,1]),9]);rm(Glioma_nor_1)
colnames(alldata)[72] <-"Glioma_nor_1"

Glioma_nor_2 <- read.table("./support_files/gene_abund_tab/SRR12518115gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Glioma_nor_2[match(alldata[,2],Glioma_nor_2[,1]),9]);rm(Glioma_nor_2)
colnames(alldata)[73] <-"Glioma_nor_2"

Glioma_nor_3 <- read.table("./support_files/gene_abund_tab/SRR12518117gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Glioma_nor_3[match(alldata[,2],Glioma_nor_3[,1]),9]);rm(Glioma_nor_3)
colnames(alldata)[74] <-"Glioma_nor_3"

Salivary_gland_nor_1 <- read.table("./support_files/gene_abund_tab/SRR13091701gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_nor_1[match(alldata[,2],Salivary_gland_nor_1[,1]),9]);rm(Salivary_gland_nor_1)
colnames(alldata)[75] <-"Salivary_gland_nor_1"

Salivary_gland_nor_2 <- read.table("./support_files/gene_abund_tab/SRR13091709gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Salivary_gland_nor_2[match(alldata[,2],Salivary_gland_nor_2[,1]),9]);rm(Salivary_gland_nor_2)
colnames(alldata)[76] <-"Salivary_gland_nor_2"

Colon_normal_1 <- read.table("./support_files/gene_abund_tab/SRR14930914gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_1[match(alldata[,2],Colon_normal_1[,1]),9]);rm(Colon_normal_1)
colnames(alldata)[77] <-"Colon_normal_1"

Colon_normal_2 <- read.table("./support_files/gene_abund_tab/SRR14930915gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_2[match(alldata[,2],Colon_normal_2[,1]),9]);rm(Colon_normal_2)
colnames(alldata)[78] <-"Colon_normal_2"

Colon_normal_3 <- read.table("./support_files/gene_abund_tab/SRR14930916gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_3[match(alldata[,2],Colon_normal_3[,1]),9]);rm(Colon_normal_3)
colnames(alldata)[79] <-"Colon_normal_3"

Colon_normal_4 <- read.table("./support_files/gene_abund_tab/SRR14930917gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_4[match(alldata[,2],Colon_normal_4[,1]),9]);rm(Colon_normal_4)
colnames(alldata)[80] <-"Colon_normal_4"

Colon_normal_5 <- read.table("./support_files/gene_abund_tab/SRR14930918gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_5[match(alldata[,2],Colon_normal_5[,1]),9]);rm(Colon_normal_5)
colnames(alldata)[81] <-"Colon_normal_5"

Colon_normal_6 <- read.table("./support_files/gene_abund_tab/SRR14930919gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_6[match(alldata[,2],Colon_normal_6[,1]),9]);rm(Colon_normal_6)
colnames(alldata)[82] <-"Colon_normal_6"

Colon_normal_7 <- read.table("./support_files/gene_abund_tab/SRR14930920gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_7[match(alldata[,2],Colon_normal_7[,1]),9]);rm(Colon_normal_7)
colnames(alldata)[83] <-"Colon_normal_7"

Colon_normal_8 <- read.table("./support_files/gene_abund_tab/SRR14930921gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_8[match(alldata[,2],Colon_normal_8[,1]),9]);rm(Colon_normal_8)
colnames(alldata)[84] <-"Colon_normal_8"

Colon_normal_9 <- read.table("./support_files/gene_abund_tab/SRR14930922gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_9[match(alldata[,2],Colon_normal_9[,1]),9]);rm(Colon_normal_9)
colnames(alldata)[85] <-"Colon_normal_9"

Colon_normal_10 <- read.table("./support_files/gene_abund_tab/SRR14930923gene_abund.txt",header = T, sep = "\t")
alldata <- cbind(alldata, Colon_normal_10[match(alldata[,2],Colon_normal_10[,1]),9]);rm(Colon_normal_10)
colnames(alldata)[86] <-"Colon_normal_10"

row.names(alldata) <- alldata[,1]; alldata <- alldata[,-c(1,2)]; alldata<-t(alldata)

alldata <- t(alldata)

alldata_nmliz <- normalize.quantiles(x=as.matrix(alldata),copy=TRUE)
colnames(alldata_nmliz) <- colnames(alldata);row.names(alldata_nmliz) <- row.names(alldata)
alldata_nmliz <- t(alldata_nmliz)

alldata_nmliz2<-log2(alldata_nmliz+1)


###
Sample_Type <- read.csv("./support_files/SampleType3.csv",header = T)  ###

gene_boxplot <- cbind(alldata, Sample_Type[match(rownames(alldata),Sample_Type[,1]),c(2,3)])
gene_boxplot1 <- cbind(alldata, Sample_Type[match(rownames(alldata),Sample_Type[,1]),c(2,3)])


##########

########

STAT3 <- gene_boxplot[,c(37,38,1)] 
STAT3_1 <- gene_boxplot1[,c(37,38,1)]
p <- ggboxplot(STAT3, x = "Type", y = "STAT3",
               color = "Group", palette = "nejm", 
               add = "none",xlab = FALSE,x.text.angle=45,size =0.9) #
p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5,method = "t.test")   

p <- ggboxplot(STAT3_1, x = "Type", y = "STAT3",
               color = "Group", palette = "nejm", 
               add = "none",xlab = FALSE,x.text.angle=45,size =0.9) #
p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5,method = "t.test")    


