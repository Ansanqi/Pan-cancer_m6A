rm(list=ls())
library(gplots)
dataraw=read.csv("mat_delNA_3log2_50.csv",header = T,row.names = 1)
data_window=read.table("./support_files/sliding_window_on_merged_tx_Ens_2020.txt",header = F,row.names = 1)

alldata <- cbind( data_window[match(rownames(dataraw), data_window$V2),c(1,3)],dataraw )
dataminus<- alldata[alldata[,2]=="-" ,]
dim(dataminus)
infominus=apply(as.matrix(dataminus), 1, function(x) unlist(strsplit(x[1], split="_"))) ; infominus=t(infominus)

newidminus=1
for (i in 2:26){
  if (infominus[i, 1]!=infominus[i-1,1]) 
    newidminus=c(newidminus, 1)
  else if (infominus[i, 1]==infominus[i-1, 1] & infominus[i, 3]!=infominus[i-1, 3])
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-1,4])>2)
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-4,4])==4 )
    newidminus=c(newidminus, newidminus[i-1]+1)
  else
    newidminus=c(newidminus, newidminus[i-1]);	
}
for (i in 27:nrow(infominus)){
  if (infominus[i, 1]!=infominus[i-1,1]) 
    newidminus=c(newidminus, 1)
  else if (infominus[i, 1]==infominus[i-1, 1] & infominus[i, 3]!=infominus[i-1, 3])
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-1,4])>=2)
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-5,4])==5 & as.numeric(infominus[i,4])-as.numeric(infominus[i-6,4])!=6 )
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-10,4])==10 & as.numeric(infominus[i,4])-as.numeric(infominus[i-11,4])!=11 )
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-15,4])==15 & as.numeric(infominus[i,4])-as.numeric(infominus[i-16,4])!=16 )
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-20,4])==20 & as.numeric(infominus[i,4])-as.numeric(infominus[i-21,4])!=21 )
    newidminus=c(newidminus, newidminus[i-1]+1)
  else if (infominus[i, 1]==infominus[i-1,1] & infominus[i, 3]==infominus[i-1,3] & as.numeric(infominus[i,4])-as.numeric(infominus[i-25,4])==25 & as.numeric(infominus[i,4])-as.numeric(infominus[i-26,4])!=26 )
    newidminus=c(newidminus, newidminus[i-1]+1)
  else
    newidminus=c(newidminus, newidminus[i-1]);	
}


newidminus2=paste(infominus[,1], newidminus, sep="_")
oldnameminus<- paste(infominus[,1],infominus[,2],infominus[,3],infominus[,4],sep = "_" )
corresponding_nameminus<- cbind(oldnameminus,newidminus2)
corresponding_name<- cbind(corresponding_nameminus,rep("-", dim(corresponding_nameminus)[1] ) )
newid2<- newidminus2

write.csv(corresponding_name,"./output_files/corresponding_name_nolung.csv")

outputminus=c()
for (i in 4:ncol(dataminus)){
  merged.maxminus=tapply(dataminus[,i], as.factor(newidminus2), function(x) max(x, na.rm=T))
  outputminus=cbind(outputminus, merged.maxminus)
}

head(outputminus)

output<- outputminus
dim(output)

ens=unlist(lapply(rownames(output), function(x) unlist(strsplit(x, split="_"))[1]))
anno=read.csv("./support_files/Ensembl_gene_type_2017_Feb4.txt", sep="\t")
genename=as.vector(anno[match(ens,anno[,1]),c(2,3)])
output2=cbind(rownames(output),ens, genename, output)
##########################################################
data2=read.csv("./support_files/out_sliding_pancancer_new.csv",row.names = 1)
data2[,-2]
newdata2<- data2[corresponding_name[,1],]
allmessage_mat<-as.matrix(cbind(corresponding_name,newdata2))
head(allmessage_mat)
dim(allmessage_mat)
allmessage_mat<-allmessage_mat[,-5]
###########################################
outputsite=c()
#for (i in 3:ncol(allmessage_mat)){
chr=tapply(allmessage_mat[,4], as.factor(newid2), function(x)unique(x))
start=tapply(allmessage_mat[,5], as.factor(newid2), function(x) min(x, na.rm=T))
end=tapply(allmessage_mat[,6], as.factor(newid2), function(x) max(x, na.rm=T))
#outputsite=cbind(output2,start,end,chr)
outputsite1=cbind(names(chr),start,end,chr)

outputsite<-cbind(  output2[match( row.names(outputsite1),output2$`rownames(output)`),],start,end,chr)


colnames(outputsite)
mat<- outputsite[,c(5:88)]##################

matcolname=read.csv("./support_files/allsample_sortcolname3batch_notissue.csv", ,header = T,na.strings = "#N/A")
colnames(mat)<- as.matrix(matcolname[1,])
rownames(mat)<-outputsite[,1]
bt<- t(mat)
colnames(mat)
newbt<-bt
a<-matrix(as.numeric(unlist(newbt)),nrow=nrow(newbt))

row.names(a)<- rownames(newbt)
output_mean<-c()
for (i in 1:ncol(a))
{
  temp<- tapply(a[,i], as.factor(rownames(a)), function(x) mean(x, na.rm=T));
  output_mean<-cbind(output_mean,temp);
  print(i)
}
colnames(output_mean)<- colnames(newbt)
output_meant<-t(output_mean)

library("preprocessCore")
library("limma")

m6amat_del<- as.data.frame(output_meant)
temp<-0
for (i in 1:nrow(m6amat_del))
{print(i);temp<- rbind(temp,colnames(sort(m6amat_del[i,],decreasing = T)[1]))}
temp <- unlist(temp)
newmat<- cbind(temp[-1,],m6amat_del,outputsite[,c(1:4,89:91)])

newmat <- apply(newmat,2,as.character)

write.csv(newmat,"./output_files/newmat_information_last_matnotissue_fixed_nolung.csv")




library(gplots)
###########sort_samples_first_merge_name
dataraw=read.csv("mat_delNA_3log2_50.csv",header = T,row.names = 1)

alldata <- cbind( data_window[match(rownames(dataraw), data_window$V2),c(1,3)],dataraw )
#################
dataplus<- alldata[alldata[,2]=="+" ,]
dim(dataplus)
dataminus<- alldata[alldata[,2]=="-" ,]
dim(dataminus)
infoplus=apply(as.matrix(dataplus), 1, function(x) unlist(strsplit(x[1], split="_")))
infoplus=t(infoplus)

for (i in 2:nrow(infoplus)){
  if (infoplus[i, 1]!=infoplus[i-1,1]) 
    newidplus=c(newidplus, 1)
  else
    newidplus=c(newidplus, newidplus[i-1]+1);	
}

newidplus2=paste(infoplus[,1], newidplus, sep="_")
oldname<- paste(infoplus[,1],infoplus[,2],infoplus[,3],infoplus[,4],sep = "_" )
corresponding_nameplus<- cbind(oldname,newidplus2)
infominus=apply(as.matrix(dataminus), 1, function(x) unlist(strsplit(x[1], split="_"))) ; infominus=t(infominus)

newidminus=1 
for (i in 2:nrow(infominus)){
  if (infominus[i, 1]!=infominus[i-1,1]) 
    newidminus=c(newidminus, 1)
  else
    newidminus=c(newidminus, newidminus[i-1]+1);	
}

newidminus2=paste(infominus[,1], newidminus, sep="_")
oldnameminus<- paste(infominus[,1],infominus[,2],infominus[,3],infominus[,4],sep = "_" )
corresponding_nameminus<- cbind(oldnameminus,newidminus2)

corresponding_name<- cbind(rbind(corresponding_nameplus ,corresponding_nameminus),c(rep("+", dim(corresponding_nameplus)[1] ),rep("-", dim(corresponding_nameminus)[1] ) ) )
newid2<- c(newidplus2,newidminus2)

write.csv(corresponding_name,"./output_files/corresponding_name_lung.csv")

dataplus<- as.matrix(dataplus)
outputplus<- dataplus[,c(4:87)];rownames(outputplus)<-newidplus2

dataminus<- as.matrix(dataminus)
outputminus<- dataminus[,c(4:87)];rownames(outputminus)<-newidminus2

head(outputminus)

output<- rbind(outputplus,outputminus)
dim(output)


ens=unlist(lapply(rownames(output), function(x) unlist(strsplit(x, split="_"))[1]))
anno=read.csv("./support_files/Ensembl_gene_type_2017_Feb4.txt", sep="\t")
genename=as.vector(anno[match(ens,anno[,1]),c(2,3)])
output2=cbind(rownames(output),ens, genename, output)
data2[,-2]
newdata2<- data2[corresponding_name[,1],]
allmessage_mat<-as.matrix(cbind(corresponding_name,newdata2))
head(allmessage_mat)
dim(allmessage_mat)
allmessage_mat<-allmessage_mat[,-5]
alldata <- cbind( data_window[match(rownames(dataraw), data_window$V2),c(1,3)],dataraw )

outputsite<-cbind(  output2[match( allmessage_mat[,2],output2$`rownames(output)`),],allmessage_mat[,5],allmessage_mat[,6],allmessage_mat[,4])
colnames(outputsite)[89] = 'start';colnames(outputsite)[90] = 'end';colnames(outputsite)[91] = 'chr'

colnames(outputsite)
mat<- outputsite[,c(5:88)]##################

matcolname=read.csv("./support_files/allsample_sortcolname3batch.csv", ,header = T,na.strings = "#N/A")
colnames(mat)<- as.matrix(matcolname[1,])
rownames(mat)<-outputsite[,1]
bt<- t(mat)
colnames(bt)<-outputsite[,1]

colnames(mat)
newbt<-bt
a<-matrix(as.numeric(unlist(newbt)),nrow=nrow(newbt))

row.names(a)<- rownames(newbt)
output_mean<-c()
for (i in 1:ncol(a))
{
  temp<- tapply(a[,i], as.factor(rownames(a)), function(x) mean(x, na.rm=T));
  output_mean<-cbind(output_mean,temp);
  print(i)
}
colnames(output_mean)<- colnames(newbt)
output_meant<-t(output_mean)

library("preprocessCore")
library("limma")

m6amat_del<- as.data.frame(output_meant)
temp<-0
for (i in 1:nrow(m6amat_del))
{print(i);temp<- rbind(temp,colnames(sort(m6amat_del[i,],decreasing = T)[1]))}
temp <- unlist(temp)
newmat<- cbind(temp[-1,],m6amat_del,outputsite[,c(1:4,89:91)])

newmat <- apply(newmat,2,as.character)
write.csv(newmat,"./output_files/newmat_information_last_matnotissue_fixed_lung.csv")
