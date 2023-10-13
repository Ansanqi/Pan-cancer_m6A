cv_ratio=function(x){  
  y=length(na.omit(x$UTR5_cv_value[x$UTR5_cv_value >=0.3 ] ) )/length(na.omit(x$UTR5_cv_value ) )
  return(y)
}

name_col=function(x){  
  y=unique(x$Group)
  return(y)
}

name_type=function(x){ 
  y=unique(x$Type)
  return(y)
}

par(mfrow=c(2,2))
barplot( cbind( c(   cv_ratio(Colon_cancer_cv),1-cv_ratio(Colon_cancer_cv)),
                c(cv_ratio(Endometrial_cancer_cv),1-cv_ratio(Endometrial_cancer_cv)),
                c(cv_ratio(Glioma_cancer_cv),1-cv_ratio(Glioma_cancer_cv)),
                c(cv_ratio(Liver_cancer_cv),1-cv_ratio(Liver_cancer_cv)),
                c(cv_ratio(Lung_cancer_cv),1-cv_ratio(Lung_cancer_cv)),
                c(cv_ratio(Ovarian_cancer_cv),1-cv_ratio(Ovarian_cancer_cv)),
                c(cv_ratio(Leukemia_cancer_cv),1-cv_ratio(Leukemia_cancer_cv)),
                c(cv_ratio(Salivary_gland_cancer_cv),1-cv_ratio(Salivary_gland_cancer_cv)),
                c(cv_ratio(Stomach_cancer_cv),1-cv_ratio(Stomach_cancer_cv)),
                c(   cv_ratio(Colon_normal_cv),1-cv_ratio(Colon_normal_cv)),
                c(cv_ratio(Endometrial_normal_cv),1-cv_ratio(Endometrial_normal_cv)),
                c(cv_ratio(Glioma_normal_cv),1-cv_ratio(Glioma_normal_cv)),
                c(cv_ratio(Liver_normal_cv),1-cv_ratio(Liver_normal_cv)),
                c(cv_ratio(Lung_normal_cv),1-cv_ratio(Lung_normal_cv)),
                c(cv_ratio(Ovarian_normal_cv),1-cv_ratio(Ovarian_normal_cv)),
                c(cv_ratio(Leukemia_normal_cv),1-cv_ratio(Leukemia_normal_cv)),
                c(cv_ratio(Salivary_gland_normal_cv),1-cv_ratio(Salivary_gland_normal_cv)),
                c(cv_ratio(Stomach_normal_cv),1-cv_ratio(Stomach_normal_cv))
) [,c(2,+2+9,8,8+9,9,9+9,3,3+9,7,7+9,1,1+9,4,4+9,5,5+9)], col = c("orange","purple") , names = cv_mat$CV_numbercol[c(2,+2+9,8,8+9,9,9+9,3,3+9,7,7+9,1,1+9,4,4+9,5,5+9)]
,las=2 ,ylab = "Fraction of peaks" )