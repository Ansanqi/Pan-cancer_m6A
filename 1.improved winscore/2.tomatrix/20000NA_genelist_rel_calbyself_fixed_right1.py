from __future__ import print_function
import pandas as pd
import os
from math import log
path = "/last_uni_peaks/"
listgene={}
filelists = os.listdir(path)
for rep in range(len(filelists)):
    bed=filelists[rep]
    for line in open(os.path.join(path, bed),'r'):
        line_1=line.rstrip('\n\r')
        temp=(line.split('\t')[3])
        listgene[temp]=1
result={}
input=open('/allmerged_winscorefile.txt','r')
for line in input:
    line=line.strip()
    a1=line.split('\t')
    result[a1[0]]=[]
    for j in range(0,84,7):          
        val=float("nan")
        if listgene.has_key(a1[0]):                 
                #ensg1=a1[j]
                #start1=a1[j+1]
		#Input_rpkm1=a1[j+2]
		#a1[j+4]=Input_rpkm_Median
		#a1[j+5]=IP_rpkm_Median
                Input_rpkm1=float("nan")
		IP_rpkm1=a1[j+3]
		if float(a1[j+2])>=5 and float(a1[j+4])>=1 and float(a1[j+5])>=1:
			Input_rpkm1=a1[j+2]
	                Input_rpkm_Median1=a1[j+4]
                	IP_rpkm_Median1=a1[j+5]
                	val =log(((float(IP_rpkm1))/(float(IP_rpkm_Median1)))/((float(Input_rpkm1))/(float(Input_rpkm_Median1)))+1,2)
		#val =(float(IP_rpkm1))/(float(Input_rpkm1))
                #winscore1=a1[j+6]
                	print(val)
        result[a1[0]].append(val)
dicfile=open('NA_first_dic_winscore_calcu_by_selfadd1_log2_fixedright.txt','w')
for key in result:
    dicfile.write(''+str(key)+'\t'+str(result[key]).strip('[').strip(']').replace(',','\t')+'\n')
    #print (result[key])
dicfile.close()
print (result)