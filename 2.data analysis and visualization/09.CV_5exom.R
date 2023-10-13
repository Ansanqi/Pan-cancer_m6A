temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])

row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

cancer_sample_name <- read.table( "./support_files/cancer_sample.txt")[,1]
normal_sample_name <- read.table( "./support_files/normal_sample.txt")[,1]

stop_conda_ID_cancer <- stop_conda_ID[,cancer_sample_name] ;stop_conda_ID_normal <- stop_conda_ID[,normal_sample_name]
UTR5_cancer <- UTR5[,cancer_sample_name]                   ;UTR5_normal <- UTR5[,normal_sample_name]
UTR3_cancer <- UTR3[,cancer_sample_name]                   ;UTR3_normal <- UTR3[,normal_sample_name]
longexon_cancer <- longexon[,cancer_sample_name]           ;longexon_normal <- longexon[,normal_sample_name]
CDS_cancer <- CDS[,cancer_sample_name]                     ;CDS_normal <- CDS[,normal_sample_name]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

stop_conda_ID_cancer_cv <- as.data.frame(apply(stop_conda_ID_cancer, 1, cal_cv)) ;stop_conda_ID_normal_cv <- as.data.frame(apply(stop_conda_ID_normal, 1, cal_cv))
stop_conda_ID_cancer_cv[,2] <- rep("Stop_codon",nrow(stop_conda_ID_cancer_cv));stop_conda_ID_cancer_cv[,3] <- rep("cancer",nrow(stop_conda_ID_cancer_cv))
stop_conda_ID_normal_cv[,2] <- rep("Stop_codon",nrow(stop_conda_ID_normal_cv));stop_conda_ID_normal_cv[,3] <- rep("normal",nrow(stop_conda_ID_normal_cv))
colnames(stop_conda_ID_cancer_cv) <- c("cv_value","Type","Group");rownames(stop_conda_ID_cancer_cv)<-NULL
colnames(stop_conda_ID_normal_cv) <- c("cv_value","Type","Group");rownames(stop_conda_ID_normal_cv)<-NULL

UTR5_cancer_cv <- as.data.frame(apply(UTR5_cancer, 1, cal_cv)) ;UTR5_normal_cv <- as.data.frame(apply(UTR5_normal, 1, cal_cv))
UTR5_cancer_cv[,2] <- rep("UTR5",nrow(UTR5_cancer_cv));UTR5_cancer_cv[,3] <- rep("cancer",nrow(UTR5_cancer_cv))
UTR5_normal_cv[,2] <- rep("UTR5",nrow(UTR5_normal_cv));UTR5_normal_cv[,3] <- rep("normal",nrow(UTR5_normal_cv))
colnames(UTR5_cancer_cv) <- c("cv_value","Type","Group");rownames(UTR5_cancer_cv)<-NULL
colnames(UTR5_normal_cv) <- c("cv_value","Type","Group");rownames(UTR5_normal_cv)<-NULL

UTR3_cancer_cv <- as.data.frame(apply(UTR3_cancer, 1, cal_cv)) ;UTR3_normal_cv <- as.data.frame(apply(UTR3_normal, 1, cal_cv))
UTR3_cancer_cv[,2] <- rep("UTR3",nrow(UTR3_cancer_cv));UTR3_cancer_cv[,3] <- rep("cancer",nrow(UTR3_cancer_cv))
UTR3_normal_cv[,2] <- rep("UTR3",nrow(UTR3_normal_cv));UTR3_normal_cv[,3] <- rep("normal",nrow(UTR3_normal_cv))
colnames(UTR3_cancer_cv) <- c("cv_value","Type","Group");rownames(UTR3_cancer_cv)<-NULL
colnames(UTR3_normal_cv) <- c("cv_value","Type","Group");rownames(UTR3_normal_cv)<-NULL

longexon_cancer_cv <- as.data.frame(apply(longexon_cancer, 1, cal_cv)) ;longexon_normal_cv <- as.data.frame(apply(longexon_normal, 1, cal_cv))
longexon_cancer_cv[,2] <- rep("Long_internal_exon",nrow(longexon_cancer_cv));longexon_cancer_cv[,3] <- rep("cancer",nrow(longexon_cancer_cv))
longexon_normal_cv[,2] <- rep("Long_internal_exon",nrow(longexon_normal_cv));longexon_normal_cv[,3] <- rep("normal",nrow(longexon_normal_cv))
colnames(longexon_cancer_cv) <- c("cv_value","Type","Group");rownames(longexon_cancer_cv)<-NULL
colnames(longexon_normal_cv) <- c("cv_value","Type","Group");rownames(longexon_normal_cv)<-NULL

CDS_cancer_cv <- as.data.frame(apply(CDS_cancer, 1, cal_cv)) ;CDS_normal_cv <- as.data.frame(apply(CDS_normal, 1, cal_cv))
CDS_cancer_cv[,2] <- rep("CDS",nrow(CDS_cancer_cv));CDS_cancer_cv[,3] <- rep("cancer",nrow(CDS_cancer_cv))
CDS_normal_cv[,2] <- rep("CDS",nrow(CDS_normal_cv));CDS_normal_cv[,3] <- rep("normal",nrow(CDS_normal_cv))
colnames(CDS_cancer_cv) <- c("cv_value","Type","Group");rownames(CDS_cancer_cv)<-NULL
colnames(CDS_normal_cv) <- c("cv_value","Type","Group");rownames(CDS_normal_cv)<-NULL

cv_alldata <- rbind(UTR5_cancer_cv,UTR5_normal_cv,CDS_cancer_cv,CDS_normal_cv,UTR3_cancer_cv,UTR3_normal_cv,longexon_cancer_cv,longexon_normal_cv,stop_conda_ID_cancer_cv,stop_conda_ID_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "cv_value",notch=TRUE,outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,1.5),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 1.5,method = "wilcox.test") + theme(text=element_text(size=32,  family="sans"))    
dev.off()




temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed(很重要).csv",header = T)

#################### 5'utr ##########
stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- UTR5[,c(1:10)]
Colon_normal <- UTR5[,c(11:20)]
Endometrial_cancer<- UTR5[,c(21:25)]
Endometrial_normal<- UTR5[,c(26:30)]
Glioma_cancer<- UTR5[,c(31:33)]
Glioma_normal<- UTR5[,c(34:36)]
Liver_cancer<- UTR5[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[3]})>0,c(37:40)]
Liver_normal<- UTR5[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[3]})>0,c(41:44)]
Lung_cancer<- UTR5[,c(45:52)]
Lung_normal<- UTR5[,c(53:56)]
Ovarian_cancer<- UTR5[,c(57:62)]
Ovarian_normal<- UTR5[,c(63:69)]
Leukemia_cancer<- UTR5[,c(70:71)]
Leukemia_normal<- UTR5[,c(72:73)]
Salivary_gland_cancer<- UTR5[,c(74:78)]
Salivary_gland_normal<- UTR5[,c(79:80)]
Stomach_cancer<- UTR5[,c(81:82)]
Stomach_normal<- UTR5[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("UTR5_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("UTR5_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer 5UTR_CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "UTR5_cv_value",outlier.size=0.1,
               color = "Group", palette = "npg",  ylim=c(0,1),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 1,method = "wilcox.test")  + theme(text=element_text(size=32,  family="sans"))   #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()

####################### CDS ###################

temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- CDS[,c(1:10)]
Colon_normal <- CDS[,c(11:20)]
Endometrial_cancer<- CDS[,c(21:25)]
Endometrial_normal<- CDS[,c(26:30)]
Glioma_cancer<- CDS[,c(31:33)]
Glioma_normal<- CDS[,c(34:36)]
Liver_cancer<- CDS[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[4]})>0,c(37:40)]
Liver_normal<- CDS[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[4]})>0,c(41:44)]
Lung_cancer<- CDS[,c(45:52)]
Lung_normal<- CDS[,c(53:56)]
Ovarian_cancer<- CDS[,c(57:62)]
Ovarian_normal<- CDS[,c(63:69)]
Leukemia_cancer<- CDS[,c(70:71)]
Leukemia_normal<- CDS[,c(72:73)]
Salivary_gland_cancer<- CDS[,c(74:78)]
Salivary_gland_normal<- CDS[,c(79:80)]
Stomach_cancer<- CDS[,c(81:82)]
Stomach_normal<- CDS[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("CDS_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("CDS_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer CDS_CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "CDS_cv_value",outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,1),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 1,method = "wilcox.test") + theme(text=element_text(size=32,  family="sans"))   #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()

###################3UTR#################


temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- UTR3[,c(1:10)]
Colon_normal <- UTR3[,c(11:20)]
Endometrial_cancer<- UTR3[,c(21:25)]
Endometrial_normal<- UTR3[,c(26:30)]
Glioma_cancer<- UTR3[,c(31:33)]
Glioma_normal<- UTR3[,c(34:36)]
Liver_cancer<- UTR3[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[4]})>0,c(37:40)]
Liver_normal<- UTR3[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[4]})>0,c(41:44)]
Lung_cancer<- UTR3[,c(45:52)]
Lung_normal<- UTR3[,c(53:56)]
Ovarian_cancer<- UTR3[,c(57:62)]
Ovarian_normal<- UTR3[,c(63:69)]
Leukemia_cancer<- UTR3[,c(70:71)]
Leukemia_normal<- UTR3[,c(72:73)]
Salivary_gland_cancer<- UTR3[,c(74:78)]
Salivary_gland_normal<- UTR3[,c(79:80)]
Stomach_cancer<- UTR3[,c(81:82)]
Stomach_normal<- UTR3[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("UTR3_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("UTR3_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer 3UTR_CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "UTR3_cv_value",outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,1),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 1,method = "wilcox.test") + theme(text=element_text(size=32,  family="sans"))    #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()

###################longexon#################


temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- longexon[,c(1:10)]
Colon_normal <- longexon[,c(11:20)]
Endometrial_cancer<- longexon[,c(21:25)]
Endometrial_normal<- longexon[,c(26:30)]
Glioma_cancer<- longexon[,c(31:33)]
Glioma_normal<- longexon[,c(34:36)]
Liver_cancer<- longexon[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[3]})>0,c(37:40)]#
Liver_normal<- longexon[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[3]})>0,c(41:44)]#
Lung_cancer<- longexon[,c(45:52)]
Lung_normal<- longexon[,c(53:56)]
Ovarian_cancer<- longexon[,c(57:62)]
Ovarian_normal<- longexon[,c(63:69)]
Leukemia_cancer<- longexon[,c(70:71)]
Leukemia_normal<- longexon[,c(72:73)]
Salivary_gland_cancer<- longexon[,c(74:78)]
Salivary_gland_normal<- longexon[,c(79:80)]
Stomach_cancer<- longexon[,c(81:82)]
Stomach_normal<- longexon[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("longexon_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("longexon_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer longexon_CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "longexon_cv_value",notch=TRUE,outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,1.5),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 1.5,method = "wilcox.test") + theme(text=element_text(size=32,  family="sans"))    #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()

###################stop_conda_ID#################


temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- stop_conda_ID[,c(1:10)]
Colon_normal <- stop_conda_ID[,c(11:20)]
Endometrial_cancer<- stop_conda_ID[,c(21:25)]
Endometrial_normal<- stop_conda_ID[,c(26:30)]
Glioma_cancer<- stop_conda_ID[,c(31:33)]
Glioma_normal<- stop_conda_ID[,c(34:36)]
Liver_cancer<- stop_conda_ID[,c(37:40)]
Liver_normal<- stop_conda_ID[,c(41:44)]
# Liver_cancer<- stop_conda_ID[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[3]})>0,c(37:40)]#
# Liver_normal<- stop_conda_ID[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[3]})>0,c(41:44)]#
Lung_cancer<- stop_conda_ID[,c(45:52)]
Lung_normal<- stop_conda_ID[,c(53:56)]
Ovarian_cancer<- stop_conda_ID[,c(57:62)]
Ovarian_normal<- stop_conda_ID[,c(63:69)]
Leukemia_cancer<- stop_conda_ID[,c(70:71)]
Leukemia_normal<- stop_conda_ID[,c(72:73)]
Salivary_gland_cancer<- stop_conda_ID[,c(74:78)]
Salivary_gland_normal<- stop_conda_ID[,c(79:80)]
Stomach_cancer<- stop_conda_ID[,c(81:82)]
Stomach_normal<- stop_conda_ID[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer Stop_codon_CV by wilcox_test_v2_不过滤肝癌.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "Stop_codon_cv_value",outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,0.8),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 0.8,method = "wilcox.test")    + theme(text=element_text(size=32,  family="sans"))   #serif  #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()


################### short #################


temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )
Shortexon<-unique( as.character(  temp_anno[temp_anno$Exon_length <400,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
Shortexon <- as.data.frame(na.omit(Shortexon));Shortexon <- cbind(Shortexon,important[match(Shortexon[,1],important[,1]),c(9:92)])


row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]
row.names(Shortexon) <- Shortexon[,1];Shortexon <- Shortexon[,-1]

Colon_cancer <- Shortexon[,c(1:10)]
Colon_normal <- Shortexon[,c(11:20)]
Endometrial_cancer<- Shortexon[,c(21:25)]
Endometrial_normal<- Shortexon[,c(26:30)]
Glioma_cancer<- Shortexon[,c(31:33)]
Glioma_normal<- Shortexon[,c(34:36)]
Liver_cancer<- Shortexon[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[4]})>0,c(37:40)]#
Liver_normal<- Shortexon[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[4]})>0,c(41:44)]#
Lung_cancer<- Shortexon[,c(45:52)]
Lung_normal<- Shortexon[,c(53:56)]
Ovarian_cancer<- Shortexon[,c(57:62)]
Ovarian_normal<- Shortexon[,c(63:69)]
Leukemia_cancer<- Shortexon[,c(70:71)]
Leukemia_normal<- Shortexon[,c(72:73)]
Salivary_gland_cancer<- Shortexon[,c(74:78)]
Salivary_gland_normal<- Shortexon[,c(79:80)]
Stomach_cancer<- Shortexon[,c(81:82)]
Stomach_normal<- Shortexon[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("Shortexon_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer shortexon_CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "Shortexon_cv_value",notch=TRUE,outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,1.5),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 1.5,method = "wilcox.test") + theme(text=element_text(size=32,  family="sans"))   #serif #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()

###############################start codon ####################################

temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

start_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==100,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

start_conda_ID <- as.data.frame(na.omit(start_conda_ID));start_conda_ID <- cbind(start_conda_ID,important[match(start_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(start_conda_ID) <- start_conda_ID[,1];start_conda_ID <- start_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- start_conda_ID[,c(1:10)]
Colon_normal <- start_conda_ID[,c(11:20)]
Endometrial_cancer<- start_conda_ID[,c(21:25)]
Endometrial_normal<- start_conda_ID[,c(26:30)]
Glioma_cancer<- start_conda_ID[,c(31:33)]
Glioma_normal<- start_conda_ID[,c(34:36)]
Liver_cancer<- start_conda_ID[,c(37:40)]
Liver_normal<- start_conda_ID[,c(41:44)]
Liver_cancer<- start_conda_ID[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[2]})>1.5,c(37:40)]#
Liver_normal<- start_conda_ID[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[2]})>1.5,c(41:44)]#
Lung_cancer<- start_conda_ID[,c(45:52)]
Lung_normal<- start_conda_ID[,c(53:56)]
Ovarian_cancer<- start_conda_ID[,c(57:62)]
Ovarian_normal<- start_conda_ID[,c(63:69)]
Leukemia_cancer<- start_conda_ID[,c(70:71)]
Leukemia_normal<- start_conda_ID[,c(72:73)]
Salivary_gland_cancer<- start_conda_ID[,c(74:78)]
Salivary_gland_normal<- start_conda_ID[,c(79:80)]
Stomach_cancer<- start_conda_ID[,c(81:82)]
Stomach_normal<- start_conda_ID[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("start_codon_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("start_codon_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

cv_alldata <- rbind(Colon_cancer_cv,Colon_normal_cv,Endometrial_cancer_cv,Endometrial_normal_cv,Glioma_cancer_cv,Glioma_normal_cv,Lung_cancer_cv,Lung_normal_cv,Ovarian_cancer_cv,Ovarian_normal_cv,Leukemia_cancer_cv,Leukemia_normal_cv,Salivary_gland_cancer_cv,Salivary_gland_normal_cv,Stomach_cancer_cv,Stomach_normal_cv,Liver_cancer_cv,Liver_normal_cv)

library(ggpubr)

pdf(file = "./output_files/nornal and cancer start_codon_CV by wilcox_test.pdf",width = 10, height = 8)
p <- ggboxplot(cv_alldata, x = "Type", y = "start_codon_cv_value",outlier.size=0.1,
               color = "Group", palette = "npg", ylim=c(0,0.8),
               add = "none",xlab = FALSE,x.text.angle=45,size =1.2) 
p + stat_compare_means(aes(group = Group),  label = "p.signif" ,label.y = 0.8,method = "wilcox.test")    + theme(text=element_text(size=32,  family="sans"))   #serif  #p + stat_compare_means(aes(group = Group),  label = "p.signif",label.y = 5) 
dev.off()

##################### fold change #####################

#start

colon_start <- mean(na.omit( Colon_cancer_cv[,1]))/mean(na.omit( Colon_normal_cv[,1]))
Endometrial_start <- mean(na.omit( Endometrial_cancer_cv[,1]))/mean(na.omit( Endometrial_normal_cv[,1]))
Glioma_start <- mean(na.omit( Glioma_cancer_cv[,1]))/mean(na.omit( Glioma_normal_cv[,1]))
Liver_start <- mean(na.omit( Liver_cancer_cv[,1]))/mean(na.omit( Liver_normal_cv[,1]))
Lung_start <- mean(na.omit( Lung_cancer_cv[,1]))/mean(na.omit(Lung_normal_cv[,1]))
Ovarian_start <- mean(na.omit( Ovarian_cancer_cv[,1]))/mean(na.omit(Ovarian_normal_cv[,1]))
Leukemia_start <- mean(na.omit(Leukemia_cancer_cv[,1]))/mean(na.omit(Leukemia_normal_cv[,1]))
Salivary_gland_start <- mean(na.omit(Salivary_gland_cancer_cv[,1]))/mean(na.omit(Salivary_gland_normal_cv[,1]))
Stomach_start <- mean(na.omit(Stomach_cancer_cv[,1]))/mean(na.omit(Stomach_normal_cv[,1]))

start_foldchange <- c(colon_start,Endometrial_start,Glioma_start,Lung_start,Ovarian_start,Leukemia_start,Salivary_gland_start,Stomach_start,Liver_start)
tissue <- c("colon","Endometrial","Glioma","Lung","Ovarian","Leukemia","Salivary_gland","Stomach","Liver")
exom <- rep("start",9)

start_foldchange <- cbind(start_foldchange,tissue,exom)
save(start_foldchange,file="./start_foldchange.RData")
rm(list = ls())


###stop
temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

stop_conda_ID<- unique( as.character(  temp_anno[temp_anno$Bin_number ==200,]$PeakID )  )
UTR5<- unique( as.character(  temp_anno[temp_anno$Bin_number <100,]$PeakID )  )
UTR3<-unique( as.character(  temp_anno[temp_anno$Bin_number >200,]$PeakID )  )
longexon<-unique( as.character(  temp_anno[temp_anno$Exon_length >400,]$PeakID )  )
CDS<- unique( as.character(  temp_anno[temp_anno$Bin_number >100 & temp_anno$Bin_number <200,]$PeakID )  )

important<- read.csv( "./support_files/newmat_information_last_matnotissue_fixed.csv",header = T)

stop_conda_ID <- as.data.frame(na.omit(stop_conda_ID));stop_conda_ID <- cbind(stop_conda_ID,important[match(stop_conda_ID[,1],important[,1]),c(9:92)])
UTR5 <- as.data.frame(na.omit(UTR5));UTR5 <- cbind(UTR5,important[match(UTR5[,1],important[,1]),c(9:92)])
UTR3 <- as.data.frame(na.omit(UTR3));UTR3 <- cbind(UTR3,important[match(UTR3[,1],important[,1]),c(9:92)])
longexon <- as.data.frame(na.omit(longexon));longexon <- cbind(longexon,important[match(longexon[,1],important[,1]),c(9:92)])
CDS <- as.data.frame(na.omit(CDS));CDS <- cbind(CDS,important[match(CDS[,1],important[,1]),c(9:92)])
row.names(stop_conda_ID) <- stop_conda_ID[,1];stop_conda_ID <- stop_conda_ID[,-1]
row.names(UTR5) <- UTR5[,1];UTR5 <- UTR5[,-1]
row.names(UTR3) <- UTR3[,1];UTR3 <- UTR3[,-1]
row.names(longexon) <- longexon[,1];longexon <- longexon[,-1]
row.names(CDS) <- CDS[,1];CDS <- CDS[,-1]

Colon_cancer <- stop_conda_ID[,c(1:10)]
Colon_normal <- stop_conda_ID[,c(11:20)]
Endometrial_cancer<- stop_conda_ID[,c(21:25)]
Endometrial_normal<- stop_conda_ID[,c(26:30)]
Glioma_cancer<- stop_conda_ID[,c(31:33)]
Glioma_normal<- stop_conda_ID[,c(34:36)]
Liver_cancer<- stop_conda_ID[,c(37:40)]
Liver_normal<- stop_conda_ID[,c(41:44)]
# Liver_cancer<- stop_conda_ID[apply(CDS[,c(37:40)] , 1, function(x){sort(x,F)[3]})>0,c(37:40)]#
# Liver_normal<- stop_conda_ID[apply(CDS[,c(41:44)] , 1, function(x){sort(x,F)[3]})>0,c(41:44)]#
Lung_cancer<- stop_conda_ID[,c(45:52)]
Lung_normal<- stop_conda_ID[,c(53:56)]
Ovarian_cancer<- stop_conda_ID[,c(57:62)]
Ovarian_normal<- stop_conda_ID[,c(63:69)]
Leukemia_cancer<- stop_conda_ID[,c(70:71)]
Leukemia_normal<- stop_conda_ID[,c(72:73)]
Salivary_gland_cancer<- stop_conda_ID[,c(74:78)]
Salivary_gland_normal<- stop_conda_ID[,c(79:80)]
Stomach_cancer<- stop_conda_ID[,c(81:82)]
Stomach_normal<- stop_conda_ID[,c(83:84)]

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}

Colon_cancer_cv <- as.data.frame(apply(Colon_cancer, 1, cal_cv));Colon_cancer_cv[,2] <- rep("Colon",nrow(Colon_cancer_cv));Colon_cancer_cv[,3] <- rep("cancer",nrow(Colon_cancer_cv));colnames(Colon_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Colon_cancer_cv)<-NULL
Colon_normal_cv <- as.data.frame(apply(Colon_normal, 1, cal_cv));Colon_normal_cv[,2] <- rep("Colon",nrow(Colon_normal_cv));Colon_normal_cv[,3] <- rep("normal",nrow(Colon_normal_cv));colnames(Colon_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Colon_normal_cv)<-NULL

Endometrial_cancer_cv <- as.data.frame(apply(Endometrial_cancer, 1, cal_cv));Endometrial_cancer_cv[,2] <- rep("Endometrial",nrow(Endometrial_cancer_cv));Endometrial_cancer_cv[,3] <- rep("cancer",nrow(Endometrial_cancer_cv));colnames(Endometrial_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Endometrial_cancer_cv)<-NULL
Endometrial_normal_cv <- as.data.frame(apply(Endometrial_normal, 1, cal_cv));Endometrial_normal_cv[,2] <- rep("Endometrial",nrow(Endometrial_normal_cv));Endometrial_normal_cv[,3] <- rep("normal",nrow(Endometrial_normal_cv));colnames(Endometrial_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Endometrial_normal_cv)<-NULL

Glioma_cancer_cv <- as.data.frame(apply(Glioma_cancer, 1, cal_cv));Glioma_cancer_cv[,2] <- rep("Glioma",nrow(Glioma_cancer_cv));Glioma_cancer_cv[,3] <- rep("cancer",nrow(Glioma_cancer_cv));colnames(Glioma_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Glioma_cancer_cv)<-NULL
Glioma_normal_cv <- as.data.frame(apply(Glioma_normal, 1, cal_cv));Glioma_normal_cv[,2] <- rep("Glioma",nrow(Glioma_normal_cv));Glioma_normal_cv[,3] <- rep("normal",nrow(Glioma_normal_cv));colnames(Glioma_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Glioma_normal_cv)<-NULL

Liver_cancer_cv <- as.data.frame(apply(Liver_cancer, 1, cal_cv));Liver_cancer_cv[,2] <- rep("Liver",nrow(Liver_cancer_cv));Liver_cancer_cv[,3] <- rep("cancer",nrow(Liver_cancer_cv));colnames(Liver_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Liver_cancer_cv)<-NULL
Liver_normal_cv <- as.data.frame(apply(Liver_normal, 1, cal_cv));Liver_normal_cv[,2] <- rep("Liver",nrow(Liver_normal_cv));Liver_normal_cv[,3] <- rep("normal",nrow(Liver_normal_cv));colnames(Liver_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Liver_normal_cv)<-NULL

Lung_cancer_cv <- as.data.frame(apply(Lung_cancer, 1, cal_cv));Lung_cancer_cv[,2] <- rep("Lung",nrow(Lung_cancer_cv));Lung_cancer_cv[,3] <- rep("cancer",nrow(Lung_cancer_cv));colnames(Lung_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Lung_cancer_cv)<-NULL
Lung_normal_cv <- as.data.frame(apply(Lung_normal, 1, cal_cv));Lung_normal_cv[,2] <- rep("Lung",nrow(Lung_normal_cv));Lung_normal_cv[,3] <- rep("normal",nrow(Lung_normal_cv));colnames(Lung_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Lung_normal_cv)<-NULL

Ovarian_cancer_cv <- as.data.frame(apply(Ovarian_cancer, 1, cal_cv));Ovarian_cancer_cv[,2] <- rep("Ovarian",nrow(Ovarian_cancer_cv));Ovarian_cancer_cv[,3] <- rep("cancer",nrow(Ovarian_cancer_cv));colnames(Ovarian_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Ovarian_cancer_cv)<-NULL
Ovarian_normal_cv <- as.data.frame(apply(Ovarian_normal, 1, cal_cv));Ovarian_normal_cv[,2] <- rep("Ovarian",nrow(Ovarian_normal_cv));Ovarian_normal_cv[,3] <- rep("normal",nrow(Ovarian_normal_cv));colnames(Ovarian_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Ovarian_normal_cv)<-NULL

Leukemia_cancer_cv <- as.data.frame(apply(Leukemia_cancer, 1, cal_cv));Leukemia_cancer_cv[,2] <- rep("Leukemia",nrow(Leukemia_cancer_cv));Leukemia_cancer_cv[,3] <- rep("cancer",nrow(Leukemia_cancer_cv));colnames(Leukemia_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Leukemia_cancer_cv)<-NULL
Leukemia_normal_cv <- as.data.frame(apply(Leukemia_normal, 1, cal_cv));Leukemia_normal_cv[,2] <- rep("Leukemia",nrow(Leukemia_normal_cv));Leukemia_normal_cv[,3] <- rep("normal",nrow(Leukemia_normal_cv));colnames(Leukemia_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Leukemia_normal_cv)<-NULL

Salivary_gland_cancer_cv <- as.data.frame(apply(Salivary_gland_cancer, 1, cal_cv));Salivary_gland_cancer_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_cancer_cv));Salivary_gland_cancer_cv[,3] <- rep("cancer",nrow(Salivary_gland_cancer_cv));colnames(Salivary_gland_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Salivary_gland_cancer_cv)<-NULL
Salivary_gland_normal_cv <- as.data.frame(apply(Salivary_gland_normal, 1, cal_cv));Salivary_gland_normal_cv[,2] <- rep("Salivary_gland",nrow(Salivary_gland_normal_cv));Salivary_gland_normal_cv[,3] <- rep("normal",nrow(Salivary_gland_normal_cv));colnames(Salivary_gland_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Salivary_gland_normal_cv)<-NULL

Stomach_cancer_cv <- as.data.frame(apply(Stomach_cancer, 1, cal_cv));Stomach_cancer_cv[,2] <- rep("Stomach",nrow(Stomach_cancer_cv));Stomach_cancer_cv[,3] <- rep("cancer",nrow(Stomach_cancer_cv));colnames(Stomach_cancer_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Stomach_cancer_cv)<-NULL
Stomach_normal_cv <- as.data.frame(apply(Stomach_normal, 1, cal_cv));Stomach_normal_cv[,2] <- rep("Stomach",nrow(Stomach_normal_cv));Stomach_normal_cv[,3] <- rep("normal",nrow(Stomach_normal_cv));colnames(Stomach_normal_cv) <- c("Stop_codon_cv_value","Type","Group");rownames(Stomach_normal_cv)<-NULL

colon_stop <- mean(na.omit( Colon_cancer_cv[,1]))/mean(na.omit( Colon_normal_cv[,1]))
Endometrial_stop <- mean(na.omit( Endometrial_cancer_cv[,1]))/mean(na.omit( Endometrial_normal_cv[,1]))
Glioma_stop <- mean(na.omit( Glioma_cancer_cv[,1]))/mean(na.omit( Glioma_normal_cv[,1]))
Liver_stop <- mean(na.omit( Liver_cancer_cv[,1]))/mean(na.omit( Liver_normal_cv[,1]))
Lung_stop <- mean(na.omit( Lung_cancer_cv[,1]))/mean(na.omit(Lung_normal_cv[,1]))
Ovarian_stop <- mean(na.omit( Ovarian_cancer_cv[,1]))/mean(na.omit(Ovarian_normal_cv[,1]))
Leukemia_stop <- mean(na.omit(Leukemia_cancer_cv[,1]))/mean(na.omit(Leukemia_normal_cv[,1]))
Salivary_gland_stop <- mean(na.omit(Salivary_gland_cancer_cv[,1]))/mean(na.omit(Salivary_gland_normal_cv[,1]))
Stomach_stop <- mean(na.omit(Stomach_cancer_cv[,1]))/mean(na.omit(Stomach_normal_cv[,1]))

stop_foldchange <- c(colon_stop,Endometrial_stop,Glioma_stop,Lung_stop,Ovarian_stop,Leukemia_stop,Salivary_gland_stop,Stomach_stop,Liver_stop)
tissue <- c("colon","Endometrial","Glioma","Lung","Ovarian","Leukemia","Salivary_gland","Stomach","Liver")
exom <- rep("stop",9)

stop_foldchange <- cbind(stop_foldchange,tissue,exom)
save(stop_foldchange,file="./stop_foldchange.RData")
rm(list = ls())

library(ggpubr)
library(ggrepel)
load(file = "./stop_foldchange.RData");load(file = "./start_foldchange.RData")
fold_change <- cbind(start_foldchange,stop_foldchange);row.names(fold_change) <- fold_change[,2]
fold_change <- as.data.frame(fold_change);fold_change <- fold_change[,-c(2,3,5,6)]
fold_change$start_foldchange <-as.numeric(fold_change$start_foldchange)
fold_change$stop_foldchange <-as.numeric(fold_change$stop_foldchange)

pdf(file = "./output_files/CV fold change.pdf",width = 5, height = 5)
p <- plot(fold_change,xlim=c(0.5,2.5),ylim=c(0.5,2.5),pch = 1,cex=1, frame = T,col ="#000080",xlab = "start codon CV-value Fold change (Cancer/Normal)", ylab = "stop codon CV-value Fold change (Cancer/Normal)")
x <- c(1,2,3,4,5);y <- c(1,2,3,4,5)
abline(lm(y ~ x), col = "grey",lwd=1.5, lty = 4)
dev.off()




colnames(fold_change)[1] <- "Fold_change"


pdf(file = "./output_files/fold change.pdf",width = 5, height = 4)
ggbarplot(fold_change, x = "tissue", y = "Fold_change", color = "exom",fill="exom",
          palette = c("orange", "purple"),x.text.angle=45,
          position = position_dodge(0.8),add = c("mean_se"))
dev.off()









