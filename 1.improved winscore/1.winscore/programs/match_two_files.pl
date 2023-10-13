#!/usr/bin/perl

#script will get name and the position-1 and their occurance from the file.sam file and put out file.stats file which looks like

#usage: occurance_sam.pl <file.sam>

die "Usage: $0  <file.sam><file2.sam><out> " unless (scalar(@ARGV)==3);

my $date =`date`;
print "$date";
print "Start to Run the Script...\n";
our %h;

@fbase=split('\.',"$ARGV[0]");

my $sam=$ARGV[0];
open (SAM, "< $sam") or die "Can not open file \"$sam\"";
while ($line=<SAM>)
 {
  chomp($line);
  $row ++;
  if ($row%1000000 eq 0) {
    print "$row\n";
  }
#Dhx9_exon_1_1	-0.969584926127501
 # if ($line=~/_range=(ch\w+):(\d+)-\d+/){
  @words = split /\t/, $line;
            $key=$words[0];#0
            $h{$key}=$words[6];#6
     #     print "$key\n" }   
}
close (SAM);
$date =`date`;
print "$date";
print "Finished File $fss was saved!\n";

$result_file=$ARGV[2];
open (RES, "> $result_file") or die "Can not open file \"$result_file\"";
my $sam2=$ARGV[1];
open (SAM2, "< $sam2") or die "Can not open file \"$sam2\"";
while ($line2=<SAM2>)
 {
  chomp($line2);
  $row ++;
  if ($row%1000000 eq 0) {
    print "$row\n";
  }
#Dhx9_exon_1_1  -0.969584926127501
 # if ($line=~/_range=(ch\w+):(\d+)-\d+/){
  @words2 = split /\t/, $line2;
  if ($h{$words2[1]}) {
  print RES "$line2\t$words2[1]\t$h{$words2[1]}\n";#1
     #     print "$key\n" 
  }   
}
close (SAM2);


