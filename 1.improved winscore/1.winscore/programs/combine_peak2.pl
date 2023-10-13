#!/usr/bin/perl

#find window within same ch and name above 2 and mark peak

#usage: combpeak.pl <input> 

die "Usage: $0  <input> " unless (scalar(@ARGV)==1);

my $date =`date`;
print "$date start to run\n";


my $fin=$ARGV[0];
my $fco=$fin.".co";
my $fpk=$fco.".pk2";
open (FIN, "< $fco") or die "Can not open file \"$fin\"";
open (RES, "> $fpk") or die "Can not open file \"$fco\"";
our %h;
our %h1;
our %h2;
our %h3;
our %h4;

for (<FIN>)
{
#chrom	txStart	txEnd	GeneSymbol	winscore	IndexWin
#chr1	155303437	155303537	Dhx9_exon_1_12	2.41231301708297	1
  chomp;
  @words=split;
  my ($chrom,$start,$end,$symbol,$score,$windex)=@words;
  if ($chrom eq "chrom")
  {print  "$_\tIndexWin\n";}
  else { 
  $h{$windex}=join("-",$score, $h{$windex});
  $h1{$windex}=join("-",$chrom, $h1{$windex});
  $h2{$windex}=join("-",$start, $h2{$windex});
  $h3{$windex}=join("-",$end, $h3{$windex});
  $h4{$windex}=join("-",$symbol, $h4{$windex});}
}
close (FIN);
$date =`date`;
print "$date finish collapsing $fin to $fco\n";

for my $key (sort keys %h)
 {
 @value=split ("\-", $h{$key});
 @value1=split ("\-", $h1{$key});
 @value2=split ("\-", $h2{$key});
 @value3=split ("\-", $h3{$key});
 @name=split ("\-", $h4{$key});
 $num=scalar(@value);
  if ($num gt 1)  {
 my $count=1;
 my $count2=1;
 for ($i=1;$i<$num;$i++)
 { if ($value[0] > $value[$i]) { $count++;}
   if ($value[$i-1] < $value[$i]) { $count2++;}
 if ($count eq $num) { print RES "$value1[0]\t$value2[0]\t$value3[0]\t$name[0]\t$value[0]\n";}
 if ($count2 eq $num) { print RES "$value1[$num-1]\t$value2[$num-1]\t$value3[$num-1]\t$name[$num-1]\t$value[$num-1]\n";}
  }
  }
  else {print RES "$value1[0]\t$value2[0]\t$value3[0]\t$name[0]\t$value[0]\n";}
  }
close(RES);


