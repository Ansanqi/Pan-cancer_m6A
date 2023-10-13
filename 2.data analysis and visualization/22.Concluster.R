
library(survival)
library(ggplot2) 
library(dplyr)
library(ggthemes)
library(ggsci)
library(rio)
library(ggstatsplot)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(pheatmap)
library(stringr)
library(factoextra)
########
alldata <- read.csv("./support_files/TCGA_31type_experssion_meanTPM5.csv",header = T,row.names = 1) 
TCGA_survival <- read.csv("./support_files/TCGA_31type_survival.csv",header = T)


alldata_d = sweep(log2(alldata+1),1, apply(log2(alldata+1),1,mean,na.rm=T))#;alldata_d <- na.omit(alldata_d)
TCGA_result_cluster <- ConsensusClusterPlus(as.matrix(alldata_d),tmyPal = c('white','white','white','#FF9D6F','#FF8040','#F75000','#F75000','#D94600','#D94600','#D94600'),maxK=15,reps=100,pItem = 0.8,pFeature = 1,clusterAlg = 'hc',distance="spearman",seed=1262118388.71279,plot="png")

TCGA_resulte_num <- as.matrix(TCGA_result_cluster[[6]]$consensusClass) ####
TCGA_resulte_num[,1] = str_pad(TCGA_resulte_num[,1], 2, side = "left", "0")
C<-c(rep("C", nrow(TCGA_resulte_num)));TCGA_resulte_num <- cbind(TCGA_resulte_num,C)
Cnew <- paste(TCGA_resulte_num[,2],TCGA_resulte_num[,1],sep = "");TCGA_resulte_num <- cbind(TCGA_resulte_num,Cnew)


TCGA_survival_mat <- b[,c(4,1,9)];rm(a);rm(b);rm(c);rm(alldata_t);rm(alldata_dist);rm(alldata_dist_ward)
TCGA_survival_mat <- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival[,3]),4])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),3])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),2])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),1])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),8])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),7])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),9])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),10])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),11])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),12])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),13])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),14])
TCGA_survival_mat<- cbind(TCGA_survival_mat,TCGA_survival[match(row.names(TCGA_survival_mat),TCGA_survival$Patient_ID),15])
colnames(TCGA_survival_mat) <- c("num","C","Cnew","Type","Sample_ID","Survival_status","Survival_month","Tumor_Stage_Code","Lymph_Node_Stage_Code","Metastasis_Stage_Code","B_cell_TIMER","T_cell_CD4_TIMER","T_cell_CD8_TIMER","Neutrophil_TIMER","Macrophage_TIMER","Myeloid_dendritic_cell_TIMER")


#display.brewer.all()
color <- brewer.pal(12,"Paired")
color1 <- brewer.pal(12,"Paired")
color2 <- brewer.pal(9,"Oranges")
color3 <- brewer.pal(9,"BuGn")
color4 <- brewer.pal(9,"Blues")
color5 <- brewer.pal(11,"RdBu")


alldata_h <- log10(alldata+1)
in_sam <-c(TCGA_survival_mat[,5])
in_sam <- na.omit(in_sam)
alldata_h <- alldata_h[,in_sam]
alldata_h <- na.omit(alldata_h)

annotation_col = data.frame(Group = TCGA_survival_mat[,3],Type = TCGA_survival_mat[,4],Tumor_Stage_Code = TCGA_survival_mat[,8] , Metastasis_Stage_Code =TCGA_survival_mat[,10],B_cell_TIMER = TCGA_survival_mat[,11],T_cell_CD4_TIMER = TCGA_survival_mat[,12],T_cell_CD8_TIMER = TCGA_survival_mat[,13],Neutrophil_TIMER = TCGA_survival_mat[,14],Macrophage_TIMER = TCGA_survival_mat[,15],Myeloid_dendritic_cell_TIMER = TCGA_survival_mat[,16]) #Lymph_Node_Stage_Code = TCGA_survival_mat[,9],#增加Type分类信息
rownames(annotation_col) = rownames(TCGA_survival_mat)
ann_colors<- list(Group = c(C01="#A6CEE3",C02="#1F78B4",C03="#B2DF8A",C04="#33A02C",C05="#FB9A99",C06="#E31A1C",C07="#FDBF6F",C08="#FF7F00",C09="#CAB2D6"),#,C10="#6A3D9A",C11="#FFFF99",C12="#B15928"
                  Tumor_Stage_Code = c(TX="#FDD0A2", T0="#FDAE6B", Tis="#FD8D3C", T1="#F16913", T2="#D94801", T3="#A63603", T4="#7F2704"),
                  Metastasis_Stage_Code = c(MX="#C6DBEF", M0="#4292C6", M1="#08306B" ),
                  Type = c(ACC="#DB7093",BLCA="#FF4500",BRCA="#FFA500",CESC="#4169E1",CHOL="#FFFF00",COAD="#0000FF",DLBC="#1E90FF",ESCA="#808000",GBM="#556B2F",HNSC="#2E8B57",KICH="#000000",KIRC="#EEE8AA",KIRP="#5F9EA0",LGG="#DCDCDC",LIHC="#7FFFAA",LUAD="#8B4513",LUSC="#BC8F8F",MESO="#FFB6C1",OV="#FFA07A",PAAD="#00FFFF",PCPG="#FF1493",PRAD="#8B0000",SARC="#BDB76B",SKCM="#008B8B",STAD="#708090",TGCT="#9370DB",THCA="#00FF00",THYM="#87CEEB",UCEC="#DC143C",UCS="#EE82EE",UVM="#800080"), #Lymph_Node_Stage_Code = c(NX="#F7FCFD", N0="#CCECE6",  N1="#66C2A4", N2="#238B45", N3="#00441B"),
                  B_cell_TIMER=c("#F7F7F7", "#FEE0B6","#FDB863","#E08214","#B35806","#7F3B08"),
                  T_cell_CD4_TIMER=c("#F7F7F7", "#D8DAEB","#B2ABD2","#8073AC","#542788","#2D004B"),
                  T_cell_CD8_TIMER=c("#F7F7F7", "#FDE0EF","#F1B6DA","#DE77AE","#C51B7D","#8E0152"),
                  Neutrophil_TIMER=c("#F7F7F7", "#D9F0D3","#A6DBA0","#5AAE61","#1B7837","#00441B"),
                  Macrophage_TIMER=c("#F7F7F7", "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"),
                  Myeloid_dendritic_cell_TIMER=c("#F7F7F7", "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061"))
bk = unique(c(seq(-2,2, length=100)))

p <- pheatmap(alldata_h,scale = "row",border = F,show_colnames = FALSE,annotation_col = annotation_col,cluster_rows = TRUE,cluster_cols = FALSE,
              show_rownames = FALSE,breaks=bk,annotation_colors = ann_colors )
dev.off()
pdf(file = "heatmap_TCGA_6con_somesam.pdf",height = 20,width =20);p;dev.off()
png(file = "heatmap_TCGA_6con.png",height = 1200,width = 1200);p;dev.off()
tiff(file = "heatmap_TCGA_6con.tiff",height = 1200,width = 1200);p;dev.off()


########

pdf(file = "survival_TCGA_6con.pdf")
plot(surv.all<- survfit(Surv(Survival_month,Survival_status == "DECEASED") ~TCGA_survival_mat$`Cnew`, data =TCGA_survival_mat),
     col=color,xlab='survival time in months',ylab='survival probabilities',lwd=2,lwd.ticks=2,box.lwd=10)
legend('topright', 
       c("C1","C2","C3","C4","C5","C6",
         signif(summary( coxph(Surv(Survival_month,Survival_status == "DECEASED") ~TCGA_survival_mat$`Cnew`, data =TCGA_survival_mat))$coefficients[1,5],3)  ),
       lty = c('solid'),col=color, bty = "n", cex = 0.8) 
dev.off()

#############
write.csv(TCGA_survival_mat,file = "./survival_mat_6con.csv")

