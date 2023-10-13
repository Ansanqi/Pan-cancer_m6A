#!/usr/bin/perl

#script will calculate the winscore 
#as winscore=log2((IP_rpkm+1)/(MedianIP_rpkm+1))/((Input_rpkm+1)/(MedianInput_rpkm+1)));

#Your input file should looks like this
#Transcript_id	GeneSymbol	Input_rpkm              IP_rpkm
#Dhx9_exon_1_1	Dhx9_exon_1_1	6.06050373581704	1.54725155153835
#Dhx9_exon_1_2	Dhx9_exon_1_2	15.2176152198618	4.64175465461505

#usage: calcuate_winscore.pl <file.txt>
die "Usage: $0  <file.txt> " unless (scalar(@ARGV)==1);
my $date =`date`;
print "$date";
print "Start to Run the Script...\n";
our %h1;
our %h2;
our %med1;
our %med2;

@fbase=split('\.',"$ARGV[0]");
my $result_file="$fbase[0]"."_winscore.txt";

my $sam=$ARGV[0];
open (SAM, "< $sam") or die "Can not open file \"$sam\"";
while ($line=<SAM>)
 {
  chomp($line);;
  $row ++;
  if ($row%1000000 eq 0) {
    print "$row\n";
  }
#Dhx9_exon_1_2  Dhx9_exon_1_2   15.2176152198618        4.64175465461505
  @words = split /\t/, $line;
  @genes = split('\_exon',"$words[0]");

   $h1{$genes[0]}=$h1{$genes[0]}."+".$words[2];    
   $h2{$genes[0]}=$h2{$genes[0]}."+".$words[3];    
}
close (SAM);
$date =`date`;
print "$date";
print "Finished File $fss was saved!\n";

our $total;
open (RES, "> $result_file") or die "Can not open file \"$result_file\"";



for my $key (sort keys %h1)
 {
# print "$key hash1 $h1{$key}\n";
# print "$key hash2 $h2{$key}\n";
 @values1 = split /\+/, $h1{$key};
 @values2 = split /\+/, $h2{$key};
# print "value1 @values1\n";
# print "value2 @values2\n";
 @array=@values1;  
 $med1{$key} = &median;
 @array=@values2;  
 $med2{$key} = &median;
# print "med1 $med1\n";
# print "med2 $med2\n";
 }

$row=0;
open (SAM, "< $sam") or die "Can not open file \"$sam\"";
while ($line=<SAM>)
 {
  chomp($line);;
  $row ++;
  if ($row eq 1) {
 @words = split /\t/, $line;
 $newwords1=$words[2]."_Median";
 $newwords2=$words[3]."_Median";
    print RES "$line\t$newwords1\t$newwords2\twinscore\n";
}
  else{
   @words = split /\t/, $line;
  @genes = split('\_exon',"$words[0]");

  $winscore=log(($words[3]+1)*($med1{$genes[0]}+1)/($med2{$genes[0]}+1)/($words[2]+1))/log(2);
  print RES "$line\t$med1{$genes[0]}\t$med2{$genes[0]}\t$winscore\n";
}
}
close (SAM);


close(RES);



sub median {
my $count = scalar (@array)+1; 
my @array = sort { $a <=> $b } @array; 
if ($count % 2) { 
 $median = $array[int($count/2)]; 
} else { 
 $median = ($array[$count/2] + $array[$count/2 - 1]) / 2; 
}
}



