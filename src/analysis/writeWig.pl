#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $bedFile $outFolder $wigP $wigN);

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("i=s"       => \$bedFile,
	    "o=s"       => \$outFolder,
	    "p=s"       => \$wigP,
	    "n=s"       => \$wigN,
	    "help"      => \$help,
             "h"        => \$help);
usage() if ($help || !$bedFile || !$outFolder);

print "writeWig.pl: started (".prettyTime().")\n";

writeWigFiles($bedFile);

print "writeWig.pl: done (".prettyTime().")\n";

sub usage {
  print STDERR "\nusage: writeWig.pl -i <file> -o <dir>\n";
  print STDERR "write wig file\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -i <file>    bed file\n";
  print STDERR " -o <file>    output folder\n";
  print STDERR " -h <file>    this (usefull) help message\n";
  print STDERR "[VERSION]\n";
  print STDERR " 12-07-2010\n";
  print STDERR "[BUGS]\n";
  print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
  print STDERR "\n";
  exit(-1);
}

sub prettyTime{
  my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
  my ($second, $minute, $hour, $dayOfMonth, $month, 
    $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  my $year = 1900 + $yearOffset;
  return "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
}

sub writeWigFiles{
  my ($file) = @_;

  open(POS, ">$outFolder/$wigP") || die "cannot open $outFolder/$wigP\n";
  open(NEG, ">$outFolder/$wigN") || die "cannot open $outFolder/$wigN\n";

  print POS "track type=wiggle_0 name=\"DARIO - read density (+)\" description=\"DARIO - read density (+)\" visibility=full\n";
  print POS "browser hide all\n";
  print POS "browser full wgRna tRNAs rnaGene refGene multiz28way\n";

  print NEG "track type=wiggle_0 name=\"DARIO - read density (-)\" description=\"DARIO - read density (-)\" visibility=full\n";
  print NEG "browser hide all\n";
  print NEG "browser full wgRna tRNAs rnaGene refGene multiz28way\n";

  my %hash = ();
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    next if(/\#/);
    my ($chr, $start, $end, $id, $expr, $strand, $freq) = split(/\s+/, $_);
    for(my $i = $start; $i <= $end; $i++){
      $hash{$chr}{$strand}{$i} += $expr;
    }
  }
  close(FILE);
  foreach my $chrom (sort {$a cmp $b} keys %hash){
    print POS "variableStep chrom=$chrom\n";
    print NEG "variableStep chrom=$chrom\n";
    foreach my $pos (sort {$a <=> $b} keys %{$hash{$chrom}{"+"}}){
      print POS "$pos\t".$hash{$chrom}{"+"}{$pos}."\n";
    }
    foreach my $pos (sort {$a <=> $b} keys %{$hash{$chrom}{"-"}}){
      print NEG "$pos\t".$hash{$chrom}{"-"}{$pos}."\n";
    }
  }
  close(POS);
  close(NEG);
}
