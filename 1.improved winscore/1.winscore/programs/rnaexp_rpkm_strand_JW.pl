#!/usr/bin/perl
#
#######################################################################################
#
# this script is used to find rna expression in certain region
#
# input files are: region_location_file, chrsize file, and rna_expression_BedGraph_files
# output file is: region_rna_expression_file
#
# RPKM normalization is applied
#
#######################################################################################


print "command line: perl rnaexp_rpkm_strand.pl annotation_file chrsize_file rna_expression_BedGraph_files output_file\n";



$raw_read_length = 50;

$location_file = $ARGV[0];
$chrsize = $ARGV[1];
for ($i=2; $i<$#ARGV; $i=$i+2) {
  push (@rnaexp_files_forward, $ARGV[$i]);
}
for ($i=3; $i<$#ARGV; $i=$i+2) {
  push (@rnaexp_files_reverse, $ARGV[$i]);
}

$outfile = $ARGV[-1];


my $totalRead="$outfile"."total_reads.txt";


# read in chromosome file
open (in, "<$chrsize");
while ($line=<in>) {
  chomp $line;
  ($chr, $size) = split /[\t+\s+]/, $line;
  push (@chrs, $chr);
}
close in;


open (in, "<$location_file");
while ($line=<in>) {
  if (!($line=~/^\#/)) {
    chomp $line;
    @data = split /\t/, $line;
    $transcript_id = $data[1];
    $strand = $data[3];

    $data[9]=~s/\,$//;
    $data[10]=~s/\,$//;

    @exon_starts = split /\,/, $data[9];
    @exon_ends = split /\,/, $data[10];
 
    $effective_length = 0;
    for ($i=0; $i<=$#exon_starts; $i++) {
      $effective_length = $exon_ends[$i] - $exon_starts[$i] + 1 + $effective_length;
    }
    
    $effective_length{$transcript_id} = $effective_length;
    $strand{$transcript_id} = $strand;

    if ( $effective_length eq 0) {
      print "$line\n";
    }
  }
}
close in;


#=head


for ($i=0; $i<=$#rnaexp_files_forward; $i++) {
  
  print "$rnaexp_files_forward[$i]\n";  
  open (in, "<$rnaexp_files_forward[$i]");
  
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    $sum{$i} = $sum{$i} + abs($data[3]*($data[2]-$data[1]));
  }
  
  close in;

  open (in, "<$rnaexp_files_reverse[$i]");
  print "$rnaexp_files_reverse[$i]\n";

  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    $sum{$i} = $sum{$i} + abs($data[3]*($data[2]-$data[1]));
  }

  close in;

}


open (out, ">$totalRead");
for ($i=0; $i<=$#rnaexp_files_forward; $i++) {
  my $actualLen=int($sum{$i}/$raw_read_length);
  print out "$i\t$actualLen\n";
}
close out;

#=cut


open (in, "<$totalRead");
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  $sum{$data[0]} = $data[1];
}
close in;



open (out, ">$outfile");
print out "Transcript_id\tGeneSymbol";
for ($i=0; $i<=$#rnaexp_files_forward; $i++) {
  $file = $rnaexp_files_forward[$i];
  $file =~s/forward//;
  print out "\t$file";
}
print out "\n";


#@chrs = ("chr1");

foreach $chr (@chrs) {

  undef %flag;
  undef %exp_forward;
  undef %exp_reverse;

  my %flag = ();
  my %exp_forward = ();
  my %exp_reverse = ();

  open (in, "<$location_file");
  while ($line=<in>) {
    if (!($line=~/^\#/)) {
      chomp $line;
      @data = split /\t/, $line;
      if ($data[2] eq $chr) {

	$data[9]=~s/\,$//;
	$data[10]=~s/\,$//;

	@exon_starts = split /\,/, $data[9];
	@exon_ends = split /\,/, $data[10];
	for ($i=0; $i<=$#exon_starts; $i++) {
	  for ($j=$exon_starts[$i]; $j<=$exon_ends[$i]; $j++) {
	    $flag{$j} = 1;
	  }
	}
      }
    }
  }
  close in;
  

  for ($i=0; $i<=$#rnaexp_files_forward; $i++) {

    $file = $rnaexp_files_forward[$i];
    print "$chr\t$file\n";
    
    open (in, "<$file");   
    $count = 0;
    while ($line=<in>) {
      chomp $line;
      @data = split /[\s+\t+]/, $line;

      if ($data[0] eq $chr) {
      
	$count ++;
	if ($count%100000 eq 0) {
	  print "$chr\t$count\t$line\n";
	}
      
	for ($j=$data[1]; $j<=$data[2]; $j++) {
	  if ($flag{$j}) {
	    $exp_forward{$j}{$i} = abs($data[3]);
	  }
	}
      }
    }
    close in;
  }
  
  for ($i=0; $i<=$#rnaexp_files_reverse; $i++) {

    $file = $rnaexp_files_reverse[$i];
    print "$chr\t$file\n";
    
    open (in, "<$file");
    $count = 0;
    while ($line=<in>) {
      chomp $line;
      @data = split /[\s+\t+]/, $line;
      
      if ($data[0] eq $chr) {
	  
	$count ++;
	if ($count%100000 eq 0) {
	  print "$chr\t$count\t$line\n";
        }
	    
	for ($j=$data[1]; $j<=$data[2]; $j++) {
	  if ($flag{$j}) {
	    $exp_reverse{$j}{$i} = abs($data[3]);
	  }
        }
      }
    }
    close in;
  }


  
  open (in, "<$location_file");
  $count = 0;
  
  while ($line=<in>) {      
    if (!($line=~/^\#/)) {
      chomp $line;
      @data = split /\t/, $line;
    
      $transcript_id = $data[1];
      $strand = $data[3];
      $gene_name = $data[12];
      
      if ($chr eq $data[2]) {
	print out "$transcript_id\t$gene_name";

	$data[9]=~s/\,$//;
	$data[10]=~s/\,$//;

	if ($strand eq "-") {

	  for ($k=0; $k<=$#rnaexp_files_reverse; $k++) {
	  
	    $rpkm = 0;
	  
	    @exon_starts = split /\,/, $data[9];
	    @exon_ends = split /\,/, $data[10];
	    for ($i=0; $i<=$#exon_starts; $i++) {
	      for ($j=$exon_starts[$i]; $j<=$exon_ends[$i]; $j++) {
		$rpkm = $rpkm + $exp_reverse{$j}{$k};
	      }
	    }

	    $rpkm = $rpkm/$effective_length{$transcript_id}*1000000000/$sum{$k}/$raw_read_length;    
	  
	    if ($rpkm) {
	      print out "\t$rpkm";
	    } 
	    else {
	      print out "\t0";
	    }
	  }
	  print out "\n";
        }
	else {
	  for ($k=0; $k<=$#rnaexp_files_forward; $k++) {

	    $rpkm = 0;

            @exon_starts = split /\,/, $data[9];
            @exon_ends = split /\,/, $data[10];
            for ($i=0; $i<=$#exon_starts; $i++) {
              for ($j=$exon_starts[$i]; $j<=$exon_ends[$i]; $j++) {
                $rpkm = $rpkm + $exp_forward{$j}{$k};
              }
            }

            $rpkm = $rpkm/$effective_length{$transcript_id}*1000000000/$sum{$k}/$raw_read_length;

            if ($rpkm) {
              print out "\t$rpkm";
            }
            else {
              print out "\t0";
            }
          }
          print out "\n";
        }
      }
    }
  }
  close in;
  
}

close out;
