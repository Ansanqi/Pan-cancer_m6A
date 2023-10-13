library(reticulate)
use_python("/python", required = T)
library(m6Aexpress)
f1 <- system.file("cancer", "SRR5194784sort_bam", package="m6Aexpress")
f2 <- system.file("cancer", "SRR5194788sort_bam", package="m6Aexpress")
f3 <- system.file("cancer", "SRR5194782_sort.bam", package="m6Aexpress")
f4 <- system.file("cancer", "SRR5194786_sort.bam", package="m6Aexpress")
f5 <- system.file("cancer", "SRR5194783sort_bam", package="m6Aexpress")
f6 <- system.file("cancer", "SRR5194787sort_bam", package="m6Aexpress")
f7 <- system.file("cancer", "SRR5194781_sort.bam", package="m6Aexpress")
f8 <- system.file("cancer", "SRR5194785_sort.bam", package="m6Aexpress")
IP_BAM <- c(f1,f2)
INPUT_BAM <- c(f5,f6)
TREATED_IP_BAM <- c(f3,f4)
TREATED_INPUT_BAM <- c(f7,f8)

##Input the gene annotation file##  
gtf <- system.file("cancer", "Homo_sapiens.GRCh38.101.gtf", package="m6Aexpress")

##Obtain the consistent peak sites##
Get_peak_infor <- Get_peakinfor(IP_BAM, INPUT_BAM,TREATED_IP_BAM, TREATED_INPUT_BAM, GENE_ANNO_GTF=gtf)

##Calling differential methylated (DM) peaks among the consistent peaks##
DM_sites_infor <- DM_detect(peak_inform=Get_peak_infor,DM_CUTOFF_TYPE="pvalue",num_ctl=2, diff_peak_pvalue=0.2)


##Calculate the methylation intensity for each gene with DM peaks##
gene_methyintensity <- gene_methy_intensity(peak_inform=DM_sites_infor,txdbinfor=NA,GENE_ANNO_GTF=gtf)


##Obtain their gene expression for INPUT samples##
get_gene_express <- Get_express_data(INPUT_BAM=c(INPUT_BAM,TREATED_INPUT_BAM ), 
                                      isPairedEnd=FALSE,GENE_ANNO_GTF = gtf,isGTFAnnotationFile=TRUE)


##Identify the differential expression gene##
obtain_DEgene <- Select_DEgene(gene_count_infor=get_gene_express,
                               cond1="control", 
                               cond2="treated",
                               num_cond1=2, 
                               num_cond2=2,
                               DIFF_GENE_cutoff_FDR=0.2,
                               DE_CUTOFF_TYPE="padj") 


##Select genes with both differential expression and differential methylation##
expr_methy_gene <- match_expr_methy(gene_count_infor=obtain_DEgene, 
                                     gene_methy_infor=gene_methyintensity,
                                    OUTPUT_DIR=NA)


##Predict m6A-reg-exp genes##
m6Aexpress_result <- m6A_Express_model(Input_file=expr_methy_gene,
                                           CUTOFF_TYPE="FDR", 
                                            FDR=0.2)


##Add differential expression and differential methylation in the result##
m6A_express_addLFC_DDM <- add_LFC_DDM(expre_methyre=m6Aexpress_result, 
                                    DE_gene=obtain_DEgene, methy_intensity=gene_methyintensity,
                                    num_cond1=2, OUTPUT_DIR=NA)
