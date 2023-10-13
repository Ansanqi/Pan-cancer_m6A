#source activate root

for i in *.bed; do echo $i; bedtools getfasta -fi /Homo_sapiens.GRCh38.dna.toplevel.fa -bed $i -fo ${i}".fa" -s -name ;done
wait
for i in *.fa; do echo $i;fasta-shuffle-letters  $i  ${i}"control.fa"  -line 1000 ;done
wait
######for i in *.fa; do echo $i;fasta-shuffle-letters  $i  ${i}"control.fa" ;done
for i in *.bed; do echo $i;findMotifs.pl   ${i}".fa"    human ${i}"analysis_output"/ -fasta   ${i}"control.fa" -rna;done
##############homer2 denovo -i drugm6aunion.fa -b drugm6aunioncontrol.fa > outputfile.txt
#source activate python
python2.7 /motif_FC_sixmodule.py HEC1a_merged.bed.fa HEC1a_merged.bed.facontrol.fa HEC1a_cancer_bedresult.txt
python2.7 /motif_FC_sixmodule.py HEPG2_merged.bed.fa HEPG2_merged.bed.facontrol.fa HEPG2_cancer_bedresult.txt
python2.7 /motif_FC_sixmodule.py iSLK_merged.bed.fa iSLK_merged.bed.facontrol.fa iSLK_cancer_bedresult.txt
python2.7 /motif_FC_sixmodule.py MOLM13_merged.bed.fa MOLM13_merged.bed.facontrol.fa MOLM13_cancer_bedresult.txt
python2.7 /motif_FC_sixmodule.py MonoMac6_merged.bed.fa MonoMac6_merged.bed.facontrol.fa MonoMac6_cancer_bedresult.txt
python2.7 /motif_FC_sixmodule.py MT4tCELL_merged.bed.fa MT4tCELL_merged.bed.facontrol.fa MT4tCELL_cancer_bedresult.txt
python2.7 /motif_FC_sixmodule.py NB4_merged.bed.fa NB4_merged.bed.facontrol.fa NB4_cancer_bedresult.txt



