#!/usr/bin/python

import os

for line in open('./inputfile/m6A_list.txt'):
    line=line.rstrip('\r\n')
    if (line[0]=='#'):
        continue
    a=line.split('\t')
    nameinput=a[0]
    output=open(nameinput+'step1.sh', 'w')
    output.write('cd /output/n')
    output.write('bedtools genomecov -bg -split -strand + -ibam /sort_bam/'+nameinput+'_sort.bam -g /programs/hg38_size_nochr.txt > ./bedgraph_strand/'+
                 nameinput+'pairedinput.plus.bedGraph\n')
    output.write('bedtools genomecov -bg -split -strand - -ibam /sort_bam/'+nameinput+'_sort.bam -g /programs/hg38_size_nochr.txt > ./bedgraph_strand/'+
                 nameinput+'pairedinput.minus.bedGraph\n')
    output.write('perl /programs/rnaexp_rpkm_strand_JW.pl /programs/sliding_window_on_merged_tx_Ens_2020_nochr.txt /programs/hg38_size_nochr.txt ./bedgraph_strand/'+
                 nameinput+'pairedinput.minus.bedGraph ./bedgraph_strand/'+nameinput+'pairedinput.plus.bedGraph '+nameinput+'input_window_rpkm.txt\n')
    nameip=a[1]
    output=open(nameip+'step1.sh', 'w')
    output.write('cd /output/n')
    output.write('bedtools genomecov -bg -split -strand + -ibam /sort_bam/'+nameip+'_sort.bam -g /programs/hg38_size_nochr.txt > ./bedgraph_strand/'+
                 nameip+'pairedip.plus.bedGraph\n')
    output.write('bedtools genomecov -bg -split -strand - -ibam /sort_bam/'+nameip+'_sort.bam -g /programs/hg38_size_nochr.txt > ./bedgraph_strand/'+
                 nameip+'pairedip.minus.bedGraph\n')
    output.write('perl /programs/rnaexp_rpkm_strand_JW.pl /programs/sliding_window_on_merged_tx_Ens_2020_nochr.txt /programs/hg38_size_nochr.txt ./bedgraph_strand/'+
                 nameip+'pairedip.minus.bedGraph ./bedgraph_strand/'+nameip+'pairedip.plus.bedGraph '+nameip+'ip_window_rpkm.txt\n')
    output.close()
    os.system('sh '+nameip+'step1.sh')
    os.system('sh '+nameinput+'step1.sh')
