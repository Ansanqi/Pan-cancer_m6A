rm(list=ls())
library('gplots')
library("multcomp")
library('WGCNA')
library('stringr')
library("dplyr")
library("ComplexHeatmap")
library(fmsb)
library(RColorBrewer)

################1motif_radar_plot#######

motif_ratio_mat_six_module <- read.table(list.files()[1])[c(2:5),c(1:4)]
tempmat_six_module<- cbind(as.character( motif_ratio_mat_six_module$V2), as.numeric( as.character( motif_ratio_mat_six_module[,4] ) ) /as.numeric( as.character( motif_ratio_mat_six_module[,3] ) )  )
tempmat_num_6module<- as.data.frame( t(tempmat_six_module)[2,] )
row.names ( tempmat_num_6module ) <- tempmat_six_module[,1]
colnames( tempmat_num_6module   ) <- as.character( str_split_fixed(c(list.files()[1]), "_bed", 2)[1] )
last_mat<- as.data.frame( t(tempmat_num_6module ) )
#############1.1getmat4module#####
get_plotmat <- function(motif_ratio_mat)
{
  tempmat<- cbind(as.character( motif_ratio_mat$V2), as.numeric( as.character( motif_ratio_mat[,4] ) ) /as.numeric( as.character( motif_ratio_mat[,3] ) )  )
  tempmat_num<- as.data.frame( t(tempmat)[2,] )
  
  row.names ( tempmat_num ) <- tempmat[,1]
  last_mat<- as.data.frame( t(tempmat_num ) )
  return(last_mat)
}
###
for (i in 2:length(list.files())) {
  motif_ratio_mattemp_6moduletemp <- read.table(list.files()[i])[c(2:5),c(1:4)]
  motif_ratio_mat_six_module<- cbind(motif_ratio_mat_six_module,motif_ratio_mattemp_6moduletemp[,c(4)])
  last_mat<- rbind(last_mat, get_plotmat(motif_ratio_mattemp_6moduletemp))
  row.names( last_mat )[i]<- as.character( str_split_fixed(c(list.files()[i]), "_bed", 2)[1] ) 
}
motif_ratio_mat_six_module<-motif_ratio_mat_six_module[,-3]
colnames( motif_ratio_mat_six_module)<- c("color","motif",row.names( last_mat))
chisq.test( matrix( as.numeric( as.matrix(motif_ratio_mat_six_module[,c(3:11)]) ),nrow = 4   )  )$p.value

maxminscale <- data.frame(GGACA=c(0.24,0),AGACT=c(0.24,0),GGACT=c(0.24,0),GAACT=c(0.24,0))
colnames(last_mat)<- colnames(maxminscale)
scaledat <- rbind( maxminscale, last_mat)


scaleb=matrix(as.numeric(as.matrix( scaledat) ),nrow=nrow(scaledat))
colnames( scaleb ) <- colnames(maxminscale)

display.brewer.all()
color <- brewer.pal(9,"Set1")

pdf("./../output_files/Ratio_m6a_4classical_motif.pdf",width = 8,height = 6 )
radarchart(as.data.frame(scaleb[c(1:11),]),axistype=1,seg=4,plty=1,cglty=4,cglcol = "grey" ,maxmin=T,pcol=color, caxislabels= c("0","6%","12%","18%","24%"),  title="Ratio of m6A motif in different tumor")

dev.off()



