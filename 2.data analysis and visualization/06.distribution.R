coding<- list.files('./forplot')
coding
colors=c("red", "brown","green","blue","purple","black","pink","yellow","orange")

#####cancer distribution######

pdf("Distribution_mostpeaks_Pancancer_cancer_revised.pdf",width = 5,height = 5)


dis_mat<-cbind(read.csv(coding[9], sep="\t", header=F)[,2],read.csv(coding[13], sep="\t", header=F)[,2],read.csv(coding[5], sep="\t", header=F)[,2],
               read.csv(coding[11], sep="\t", header=F)[,2],read.csv(coding[18], sep="\t", header=F)[,2],read.csv(coding[14], sep="\t", header=F)[,2],
               read.csv(coding[15], sep="\t", header=F)[,2],read.csv(coding[7], sep="\t", header=F)[,2],read.csv(coding[1], sep="\t", header=F)[,2])
dis_mat <- dis_mat[-30,]##outliet
dis_matsum<- apply(dis_mat , 2, sum)
library(reshape2)
dis_mat. <- dis_mat; colnames(dis_mat.) <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9")
dis_mat. <- melt(dis_mat,id.vars=row.names(dis_mat))
AOV. <- aov(value~Var2, dis_mat.);summary(AOV.)[[1]][1,5]

plot(1:30, rep(0,30), type="n", ylim=c(0, .2),  ylab="% m6A peaks per transcript",main = "Distribution of m6A", xlab="")
for (i in 1:9){
  
  lines(dis_mat[,i]/dis_matsum[i], col=colors[i] )
}
text(c(5, 15, 25), -0.045, labels=c("5'UTR", "CDS", "3'UTR"), adj=0.5, xpd=T)
abline(v=10, lty=3)
abline(v=20, lty=3)
legend("topleft", legend=c("Colon cancer","Endometrial cancer", "Glioma","Leukemia","Liver cancer","Lung cancer","Ovarian cancer","Salivary cancer","Stomach cancer"), col=colors, fil=colors, border=colors,bty="n")
dev.off()

#####normal distribution######

pdf("Distribution_mostpeaks_Pancancer_normal_revised.pdf",width = 5,height = 5)
dis_mat<-cbind(read.csv(coding[10], sep="\t", header=F)[,2],read.csv(coding[12], sep="\t", header=F)[,2],read.csv(coding[6], sep="\t", header=F)[,2],
               read.csv(coding[4], sep="\t", header=F)[,2],read.csv(coding[17], sep="\t", header=F)[,2],read.csv(coding[2], sep="\t", header=F)[,2],
               read.csv(coding[16], sep="\t", header=F)[,2],read.csv(coding[8], sep="\t", header=F)[,2],read.csv(coding[3], sep="\t", header=F)[,2]) 
dis_mat <- dis_mat[-30,]##outliet######################
dis_matsum<- apply(dis_mat , 2, sum)

plot(1:30, rep(0,30), type="n", ylim=c(0, .2),  ylab="% m6A peaks per transcript",main = "Distribution of m6A", xlab="")
for (i in 1:9){
  
  lines(dis_mat[,i]/dis_matsum[i], col=colors[i] )
}
text(c(5, 15, 25), -0.045, labels=c("5'UTR", "CDS", "3'UTR"), adj=0.5, xpd=T)
abline(v=10, lty=3)
abline(v=20, lty=3)
legend("topleft", legend=c("Colon normal","Endometrial normal", "Glioma normal","PBMC","Liver normal","Lung normal","Ovarian normal","Salivary normal","Stomach normal"), col=colors, fil=colors, border=colors,bty="n")
dev.off()
