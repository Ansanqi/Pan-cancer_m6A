#coding=utf-8 
from __future__ import print_function
import pandas as pd
import os
a=list()
path = "/last_uni_peaks/"
listgene={}
filelists = os.listdir(path)
for rep in range(len(filelists)):
    bed=filelists[rep]
    for line in open(os.path.join(path, bed),'r'):
        line_1=line.rstrip('\n\r')
        temp=(line.split('\t')[3])
        listgene[temp]=1
#listgene=list(set(a))
result={}
input=open('/def_winscore/NA_first_dic_winscore_calcu_by_selfadd1_log2_fixedright.txt','r')
for line in input:
        line=line.strip()
        a1=line.split('\t')            
        if listgene.has_key(a1[0]):                 
                #ensg1=a1[j]
                #start1=a1[j+1]
                #Input_rpkm1=a1[j+2]
                #IP_rpkm1=a1[j+3]
                #Input_rpkm_Median1=a1[j+4]
                #IP_rpkm_Median1=a1[j+5]
                #val =(float(IP_rpkm1) + 1)/(float(Input_rpkm1) + 1)
                #winscore1=a1[j+6]
                #print(winscore1)
                #print (winscore1)
                #break
            result[a1[0]]=a1
            #print(a1)
dicfile=open('NA_second_dic_winscore_byself_caculate2add1_log2_fixed_right.txt','w')
for key in result:
    #innerkey,result[key][innerkey]
    dicfile.write(''+str(key)+'\t'+str(result[key]).strip('[').strip(']').replace(',','\t')+'\n')
    #print (result[key])
dicfile.close()
#print 'ins: %s'%(s1.intersection(s2))  
#print 'uni: %s'%(s1.union(s2))  
#print 'dif: %s'%(s1.difference(s2).union(s2.difference(s1))) 
print (result)
