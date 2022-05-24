#!/usr/bin/perl 

#!/usr/local/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $clusters $annotations $overlap $sizeCutoff $printOpt);


# -----------------------------------------------------------------------------
# OPTIONS

$overlap = -1;
$sizeCutoff = -1;
$printOpt = 0;

GetOptions ("c=s"       => \$clusters,
	    "a=s"       => \$annotations,
	    "o=s"       => \$overlap,
	    "s=s"       => \$sizeCutoff,
	    "p=s"       => \$printOpt,
	    "help"      => \$help,
             "h"        => \$help);

usage() if ($help || !$clusters || !$annotations);


my %ncRNAs = ();

my $c = 0;
open(NCRNAS, "<$annotations") || die "cannot open $annotations\n";
while(<NCRNAS>){
  chomp;
  my @line = split(/\s+/, $_);
  $ncRNAs{$c}{chrom} = $line[0];
  $ncRNAs{$c}{start} = $line[1];
  $ncRNAs{$c}{end} = $line[2];
  $ncRNAs{$c}{strand} = $line[5];
  $ncRNAs{$c}{id} = $line[3];
  $ncRNAs{$c}{source} = $line[7];
  $ncRNAs{$c}{type} = $line[6];
  $c++;
}
close(NCRNAS);

my $s = 0;
open(FILE, "<$clusters") || die "cannot open $clusters\n";
while(<FILE>){
  chomp;
  next if(/\#/);
  if(/^\>/){
    my @line = split(/\s+/, $_);
    my $flag = "a"; $s = 0;
    foreach my $c (keys %ncRNAs){
      next if($ncRNAs{$c}{chrom} ne $line[1]);
      next if($ncRNAs{$c}{strand} ne $line[4]);

     if($sizeCutoff != -1){ next if($sizeCutoff * ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1) < ($line[3] - $line[2] + 1)); }
 
     if($line[2] >= $ncRNAs{$c}{start} && $line[2] <= $ncRNAs{$c}{end} ||
	 $line[3] >= $ncRNAs{$c}{start} && $line[3] <= $ncRNAs{$c}{end} ||
	 $line[2] >= $ncRNAs{$c}{start} && $line[3] <= $ncRNAs{$c}{end} ||
	 $line[2] <= $ncRNAs{$c}{start} && $line[3] >= $ncRNAs{$c}{end}){

	my $sum = 0;
	for(my $x = $ncRNAs{$c}{start}; $x <= $ncRNAs{$c}{end}; $x++){
	  if($x >= $line[2] && $x <= $line[3]){$sum++;}
	}
	if($overlap != -1){next if($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1) < $overlap);}

	$flag = "$_\t$ncRNAs{$c}{id}\t$ncRNAs{$c}{type}";
	$s = 1;
      }
    }
    if($flag eq "a" && $printOpt == 0){
      print "$_\tn\/a\tn\/a\n";
    }
    elsif($flag ne "a"){
      print "$flag\n";
    }
  }
  elsif($s == 1){
    print "$_\n";
  }
  elsif($printOpt == 0){
    print "$_\n";
  }
}
close(FILE);


sub usage {
  print STDERR "\nusage: flagKnownClusters.pl -c <file> -a <file> [OPTIONS]\n";
  print STDERR "flag clusters overlapping with annotated ncRNAs\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -c <file>    cluster file\n";
  print STDERR " -a <file>    annotation file\n";
  print STDERR "[OPTIONS]\n";
  print STDERR " -o <file>    min overlap in % (default = off)\n";
  print STDERR " -s <file>    size cutoff, i.e. 1.5 times of the length (default = off)\n";
  print STDERR " -p <int>     print option (0 = all; 1 = only clusters with annotation)\n";
  print STDERR "[VERSION]\n";
  print STDERR " 11-24-2010\n";
  print STDERR "[BUGS]\n";
  print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
  print STDERR "\n";
  exit(-1);
}



