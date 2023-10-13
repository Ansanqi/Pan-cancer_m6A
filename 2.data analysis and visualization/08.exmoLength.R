temp_anno<- read.table( "./support_files/Bin_annotation_of_all_merged_peaks.txt",sep = "\t",header = T)
getwd()

cancerpeakbed <-read.table("./support_files/cancer_allpeak_newname.txt")
cancer_annomat<-temp_anno[match( cancerpeakbed$V4 ,temp_anno$PeakID),]

normalpeakbed <-read.table("./support_files/normal_allpeak_newname.txt")
normal_annomat<-temp_anno[match( normalpeakbed$V4 ,temp_anno$PeakID),]

cancer_length <- as.numeric(na.omit( cancer_annomat$Exon_length )) 
normal_length <- as.numeric(na.omit( normal_annomat$Exon_length )) 
p_value <- wilcox.test(cancer_length,normal_length)[["p.value"]]

pdf("./output_file/exonlengthdistributioncorlor.pdf",width = 5,height = 5 )
plot(density(log2(as.vector( na.omit( cancer_annomat$Exon_length ) ))),ylim=c(0,0.4)  ,col="#FF4500",main = "Density of target exon longth" ,xlab = "Log 2 internal exon length(bp)") 
lines(density(log2(as.vector( na.omit( normal_annomat$Exon_length ) ))),col="#008080")
legend("topleft", legend=c("Cancer","Normal"), col=c("#FF4500","#008080"), fil=c("#FF4500","#008080", border=c("#FF4500","#008080",bty="n")))
legend("topright", legend=paste("wilcox.test p=",p_value) )
dev.off()



