#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $bedFile $outFolder $ncRNAFile $outFile $wigP $wigN $RNAzFile $species);
my $absReadCnt = 0; my %RNAz = ();

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("b=s"       => \$bedFile,
	    "o=s"       => \$outFolder,
	    "a=s"       => \$ncRNAFile,
	    "f=s"       => \$outFile,
	    "p=s"       => \$wigP,
	    "n=s"       => \$wigN,
	    "r=s"       => \$RNAzFile,
	    "s=s"       => \$species,
	    "help"      => \$help,
             "h"        => \$help);
usage() if ($help || !$bedFile || !$ncRNAFile || !$outFolder || !$outFile || !$wigP || !$wigN || !$species);

print "getExpression.pl: started (".prettyTime().")\n";


getReadCnt();
my %ncRNAs = getncRNAs($ncRNAFile);
my %reads  = getReads($bedFile);
getExpression();
if($RNAzFile){%RNAz = readRNAz($RNAzFile);}
printExpression(%ncRNAs);

print "getExpression.pl: done (".prettyTime().")\n";

sub usage {
  print STDERR "\nusage: getExpression.pl -i <file> -o <dir>\n";
  print STDERR "check bed file\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -b <file>    bed file\n";
  print STDERR " -a <file>    ncRNA file\n";
  print STDERR " -o <file>    output folder\n";
  print STDERR " -f <file>    output file\n";
  print STDERR " -p <file>    wig file (+)\n";
  print STDERR " -n <file>    wig file (-)\n";
  print STDERR " -r <file>    RNAz file (-)\n";
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

sub getncRNAs{
  my ($file) = @_;
  my %hash = (); my $c = 0;
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    next if(/\#/);
    my ($chr, $start, $end, $id, $score, $strand, $type, $source) = split(/\s+/, $_);  
    next if(!$id);
    $hash{$c}{line} = $_;
    $hash{$c}{chr} = $chr;
    $hash{$c}{start} = $start;
    $hash{$c}{end} = $end;
    $hash{$c}{strand} = $strand;
    $hash{$c}{id} = $id;
    $hash{$c}{source} = $source;
    $hash{$c}{type} = $type;
    $hash{$c}{score} = $score;
    $hash{$c}{expr} = 0;
    $hash{$c}{reads} = 0;
    $hash{$c}{tags} = 0;
    $c++;
  }
  close(FILE);
  return %hash;
}

sub getReads{
  my ($file) = @_;
  my %hash = ();
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    next if(/\#/);
    my ($chrom, $start, $end, $id, $expr, $strand, $reads) = split(/\s+/,$_);
    $hash{$chrom}{$strand}{$start}{$end}{expr} = $expr;
    $hash{$chrom}{$strand}{$start}{$end}{reads} = $reads;
  }
  close(FILE);
  return %hash;
}

sub getExpression{
  foreach my $c (keys %ncRNAs){
    my $readCnt = 0;
    my $expr = 0;
    my $normReadCnt = 0;
    foreach my $chrom (keys %reads){
      next if($chrom ne $ncRNAs{$c}{chr});
      foreach my $strand (keys %{$reads{$chrom}}){
	next if($strand ne $ncRNAs{$c}{strand});
	foreach my $start (sort {$a <=> $b} keys %{$reads{$chrom}{$strand}}){
	  next if($start > $ncRNAs{$c}{end});
	  foreach my $end (keys %{$reads{$chrom}{$strand}{$start}}){
	    next if($end < $ncRNAs{$c}{start});
	    for(my $i = $start; $i <= $end; $i++){
	      if($i >= $ncRNAs{$c}{start} && $i <= $ncRNAs{$c}{end}){
		$expr += $reads{$chrom}{$strand}{$start}{$end}{expr} / ($end - $start + 1);
	      }
	    }
	    $readCnt += $reads{$chrom}{$strand}{$start}{$end}{reads};
	    $normReadCnt += $reads{$chrom}{$strand}{$start}{$end}{expr};
	  }
	}
      }
    }
    $ncRNAs{$c}{expr} = $expr;
    $ncRNAs{$c}{normReadCnt} = $normReadCnt;
    $ncRNAs{$c}{readCnt} = $readCnt;
  }
}

sub readRNAz{
  my ($file) = @_;
  my %hash = ();
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    next if(/\#/);
    my ($chrom, $start, $end, $id, $score, $strand) = split(/\s+/,$_);
    for(my $i = $start; $i <= $end; $i+=10){
      $hash{$chrom}{$i} = 1;
    }
    $hash{$chrom}{$end} = 1;
  }
  close(FILE);
  return %hash;
}

sub printExpression{
  my (%ncRNAs) = @_;
  open(OUT, ">$outFile") || die "cannot open $outFile";
  foreach my $c (keys %ncRNAs){
    next if($ncRNAs{$c}{expr} == 0);
    my $UCSC = "";
    my $folderID = $outFolder;
    $folderID =~ s/\/scratch\/dario\/computations\///; $folderID =~ s/\///;

    # Get the link to show expression
    if($species =~ /^ath/) { # create specific link for arabidopsis
        if($ncRNAs{$c}{strand} eq "+"){
            $UCSC = "http://chualab.rockefeller.edu/cgi-bin/gb2/gbrowse/arabidopsis/?start=".($ncRNAs{$c}{start} - 50).";stop=".($ncRNAs{$c}{end} + 50).";ref=$ncRNAs{$c}{chr};eurl=http://dario.bioinf.uni-leipzig.de/result/$folderID/$wigP";
        }else{
            $UCSC = "http://chualab.rockefeller.edu/cgi-bin/gb2/gbrowse/arabidopsis/?start=".($ncRNAs{$c}{start} - 50).";stop=".($ncRNAs{$c}{end} + 50).";ref=$ncRNAs{$c}{chr};eurl=http://dario.bioinf.uni-leipzig.de/result/$folderID/$wigN";
        }
    }else { # normal UCSC string
        if($ncRNAs{$c}{strand} eq "+"){
            $UCSC = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=$species&position=$ncRNAs{$c}{chr}:".($ncRNAs{$c}{start} - 50)."-".($ncRNAs{$c}{end} + 50)."&hgct_customText=http://dario.bioinf.uni-leipzig.de/result/$folderID/$wigP";
        }else{
            $UCSC = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=$species&position=$ncRNAs{$c}{chr}:".($ncRNAs{$c}{start} - 50)."-".($ncRNAs{$c}{end} + 50)."&hgct_customText=http://dario.bioinf.uni-leipzig.de/result/$folderID/$wigN";
        }
    }


    my $expr = ($ncRNAs{$c}{expr} / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1)) / $absReadCnt * 1000000;
    $expr = sprintf "%.2e", $expr;

    my $readCnt = $ncRNAs{$c}{readCnt};
    1 while $readCnt =~ s/^(-?\d+)(\d{3})/$1,$2/;

    my $normReads = $ncRNAs{$c}{normReadCnt};
    $normReads = sprintf "%.2f", $normReads;
    1 while $normReads =~ s/^(-?\d+)(\d{3})/$1,$2/;

    my $RNAzFlag = 0;
    for(my $i = $ncRNAs{$c}{start}; $i <= $ncRNAs{$c}{end}; $i++){
      if(exists($RNAz{$ncRNAs{$c}{chr}}{$i})){ 
	$RNAzFlag = 1; last;
      }
    }
    
    print OUT "$ncRNAs{$c}{chr}\t$ncRNAs{$c}{start}\t$ncRNAs{$c}{end}\t$ncRNAs{$c}{id}\t$ncRNAs{$c}{score}\t$ncRNAs{$c}{strand}\t$ncRNAs{$c}{type}\t$expr\t$readCnt\t$normReads\t$UCSC\t$RNAzFlag\n";
  }
  close(OUT);
}

sub getReadCnt{
  open(FILE, "<$outFolder\/upload.info") || die "cannot open $outFolder\/upload.info\n";
  while(<FILE>){
    chomp;
    next if(/\#/);
    my ($type,$count) = split(/\:/,$_);
    if($type eq "reads"){
      $absReadCnt = $count;
    }
  }
  close(FILE);
}
