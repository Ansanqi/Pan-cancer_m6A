#!/usr/bin/python
import sys

strand_info={}
for line in open('sliding_window_on_merged_tx_Ens_2020_nochr.txt'):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    strand_info[a[1]]=a[3]

output=open(sys.argv[2], 'w')
for line in open(sys.argv[1]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    strand=strand_info[a[3]]
    output.write(line+'\t'+strand+'\n')
    
