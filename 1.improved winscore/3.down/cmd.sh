for i in *.bedchr.bed; do echo $i;findMotifsGenome.pl $i hg38 ${i}"analysis_output/" -size 200 -mask -rna;done
