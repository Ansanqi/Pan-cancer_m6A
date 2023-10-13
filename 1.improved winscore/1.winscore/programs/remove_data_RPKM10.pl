#!/opt/ActivePerl/bin/perl
die "Usage: $0  <file.sorted> <S1.file> <res.out> " unless (scalar(@ARGV)==2);
my $sorted=$ARGV[0];
my $res=$ARGV[1];


open (INPUT, "< $sorted") or die "Can not open file \"$sorted100\"";

open (RES, "> $res") or die "Can not open file \"$res\"";
my $count=0;
foreach $line (<INPUT>)
{
  chomp($line);
  
  @values=split('\t', $line);
#  print "values2$values[2]\tvalues3$values[3]\n";
  if (( $values[2] > 6 ) or ($values[3] > 6))
  {print RES "$line\n";
  $count++;
   }
}
print "$count\n";
close(INPUT);
close (RES);

