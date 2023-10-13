setwd("/def_winscore")
mat<- read.table("NA_second_dic_winscore_byself_caculate2add1_log2_fixed_right.txt",sep="\t",header = F,stringsAsFactor = FALSE,row.names = 1 )
#second_count_cause_
mat1<- mat[2,]
for (i in c(3:dim(mat)[1])){if(sum(is.na(mat[i,]))<= 10){ temp <- mat[i,];mat1 <- rbind(mat1,temp);}}
write.csv(mat1,"mat_delNA_3log2.csv")