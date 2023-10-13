#!/usr/bin/perl

#find window within same ch and name above 2 and mark peak

#usage: combpeak.pl <input> 

die "Usage: $0  <input> " unless (scalar(@ARGV)==1);

my $date =`date`;
print "$date start to run\n";


my $fin=$ARGV[0];
my $fco=$fin.".co";
my $fpk=$fco.".pk";
open (FIN, "< $fin") or die "Can not open file \"$fin\"";
open (FCO, "> $fco") or die "Can not open file \"$fco\"";


my ($ch1,$ch2,$ch3,$st1,$st2,$st3,$en1,$en2,$en3,$sy1,$sy2,$sy3,$sc1,$sc2,$sc3,$count,$status,$inw,$in1,$in2,$in3,$pk1,$pk2,$pk3);
$inw = 1;
for (<FIN>)
{

  chomp;
  @words=split;
  my ($chrom,$start,$end,$symbol,$score)=@words;
  if ($chrom eq "chrom")
  {print FCO "$_\tIndexWin\n";}
  elsif ($score >2)
  {
   print FCO "$_\t$inw\n";
   $status = "big";
  }
  elsif ($status eq "big")
  {
   $status = "small";
   $inw++;
  }
  #print "$chrom,$start,$end,$symbol,$score\n";
  

}
close (FCO);
close (FIN);
$date =`date`;
print "$date finish collapsing $fin to $fco\n";

open (FCO, "< $fco") or die "Can not open file \"$fco\"";
open (FPK, "> $fpk") or die "Can not open file \"$fpk\"";
for (<FCO>)
{
  $count++;
  chomp;
  @words=split;
  ($ch3,$st3,$en3,$sy3,$sc3,$in3)=@words;
  if ($ch3 eq "chrom")
  {print FPK "$ch3\t$st3\t$en3\t$sy3\t$sc3\tPeak\n";$count--;
   #print  "$ch3\t$st3\t$en3\t$sy3\t$sc3\tPeak\n";
   next;}
  if (($ch1 eq $ch2) and ($ch2 eq $ch3) and ($in1 == $in2) and ($in2 == $in3))
  {
    	#chrom	txStart	        txEnd	GeneSymbol	winscore
    	#chr1	155302887	155302987	Dhx9_exon_1_1	-0.969584926127501
	#chr1	155302937	155303037	Dhx9_exon_1_2	-1.02210110306143
	#chr1	155302987	155303087	Dhx9_exon_1_3	-1.11310493566165
	#chr1	155303037	155303137	Dhx9_exon_1_4	0.544260897707476
    my ($sym1,$sym2,$sym3)=($sy1,$sy2,$sy3);
    $sym1 =~ s/(.*_\d+)_\d+$/$1/;
    $sym2 =~ s/(.*_\d+)_\d+$/$1/;
    $sym3 =~ s/(.*_\d+)_\d+$/$1/;
    #print "$sym1,$sym2,$sym3\n";
    if (($sym1 eq $sym2) and ($sym3 eq $sym2))
    {
        if ($sc2>$sc1 and $sc2>$sc3)
        {$pk2="Y";} 
    }
  }
  if ($count == 3)
  {print FPK "$ch1\t$st1\t$en1\t$sy1\t$sc1\t$pk1\n";$count--;}
  ($ch1,$st1,$en1,$sy1,$sc1,$in1,$pk1)=($ch2,$st2,$en2,$sy2,$sc2,$in2,$pk2);
  ($ch2,$st2,$en2,$sy2,$sc2,$in2,$pk2)=($ch3,$st3,$en3,$sy3,$sc3,$in3,$pk3);
  

}
while ($count >0)
{
  print FPK "$ch1\t$st1\t$en1\t$sy1\t$sc1\t$pk1\n";
  ($ch1,$st1,$en1,$sy1,$sc1,$in1,$pk1)=($ch2,$st2,$en2,$sy2,$sc2,$in2,$pk2);
  $count--;
}
close (FCO);
close (FPK);
$date =`date`;
print "$date finish picking peak $fco to $fpk\n";
