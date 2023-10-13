#!/usr/bin/python

import os

for line in open('./inputfile/m6A_list.txt'):
    line=line.rstrip('\r\n')
    #if (line[0]=='#'):
    #    continue
    a=line.split('\t')
    nameinput=a[0]
    nameip=a[1]
    
    output=open(nameip+nameinput+'_step2.sh', 'w')
    output.write('cd /output\n')
    output.write('python /programs/merge_input_eluate_win_RPKM.py /output/'+nameinput+'input_window_rpkm.txt /output/'+nameip+'ip_window_rpkm.txt '+nameinput+'_input'+nameip+'ip_and_eluate_win_rpkm.txt\n')
    output.write('perl /programs/calculate_winscore.pl '+nameinput+'_input'+nameip+'ip_and_eluate_win_rpkm.txt\n')
    output.write('perl /programs/remove_data_RPKM10.pl '+nameinput+'_input'+nameip+'ip_and_eluate_win_rpkm_winscore.txt '+nameinput+'_input'+nameip+'ip_winscore_filtered.txt\n')
    output.write('perl /programs/match_two_files.pl '+nameinput+'_input'+nameip+'ip_winscore_filtered.txt /programs/sliding_window_on_merged_tx_Ens_2020_nochr.txt '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched.txt\n')
    output.write('cut -f 3,5,6,17,18 '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched.txt > '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated.txt\n')
    output.write('perl /programs/combine_peak.pl '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated.txt\n')
    output.write('perl /programs/combine_peak2.pl '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated.txt\n')
    output.write('grep -P "\tY$" '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated.txt.co.pk | cut -f 1-5 > '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated_greped.txt\n')
    output.write('cat '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated_greped.txt '+nameinput+'_input'+nameip+'ip_winscore_filtered_matched_formated.txt.co.pk2  > '+nameinput+'_input'+nameip+'ip_peaklist.txt\n')
    output.write('python /programs/add_strand_to_final_peaks.py '+nameinput+'_input'+nameip+'ip_peaklist.txt '+nameinput+'_input'+nameip+'ipfinal_peaks_with_strand.bed')
    output.close()
    os.system('sh '+nameip+nameinput+'_step2.sh')
