import pandas as pd
import numpy as np
import os
pathDir =  os.listdir('/all_batchpeak_list/Pan_cancer_m6A/')
i = 0
for i in range(0,len(pathDir)):
    child = os.path.join('%s%s' % ('/all_batchpeak_list/Pan_cancer_m6A/', pathDir[i]))
    print child
    data = pd.read_table(child,header=None,encoding='gb2312',delim_whitespace=True)
    data_sorted=data.sort_values(map(int,str(012)), ascending=False)
    IsDuplicated = data_sorted.loc[:,map(int,str(012))].duplicated()
    print IsDuplicated 
    datanew = data_sorted.drop_duplicates(map(int,str(012)))
#    print datanew[0]
    datanew.to_csv(pathDir[i],header=False,encoding='gb2312',sep ='\t',index=False)
#    f = open(pathDir[i],'w')
#    f.writelines(str(datanew))
#    f.close()
