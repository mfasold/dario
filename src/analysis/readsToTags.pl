#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $bedFile $outFile $summaryFile $lengthDistrFile $multiMapFile %unique %tags %entries);

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("i=s"       => \$bedFile,
            "o=s"       => \$outFile,
            "s=s"       => \$summaryFile,
            "l=s"       => \$lengthDistrFile,
            "m=s"       => \$multiMapFile,
            "help"      => \$help,
            "h"         => \$help);
usage() if ($help || !$bedFile || !$outFile || !$summaryFile || !$lengthDistrFile || !$multiMapFile);

# -----------------------------------------------------------------------------
# MAIN

print "readsToTags.pl: started (".prettyTime().")\n";

my $result = check($bedFile);
writeNormalizedBed($bedFile, $outFile);

print "checkUpload.pl: done (".prettyTime().")\n";

# -----------------------------------------------------------------------------
# FUNCTIONS

sub usage {
  print STDERR "\nusage: readsToTags.pl -i <file> -s <file> -l <file> -m <file> -o <file>\n";
  print STDERR "merge reads to tags\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -i <file>    bed file with reads\n";
  print STDERR " -s <file>    summary file\n";
  print STDERR " -l <file>    length distribution file\n";
  print STDERR " -m <file>    multiple mapping distribution file\n";
  print STDERR " -o <file>    bed file with tags\n";
  print STDERR " -h <file>    this (usefull) help message\n";
  print STDERR "[VERSION]\n";
  print STDERR " 06-08-2010\n";
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

sub check{
  my ($file) = @_;
  my $line = 0;
  my %length = ();
  my $readCount = 0;
  my $tagCount = 0;
  my $entryCount = 0;

  # FIND MULTIPLE MAPPINGS
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    next if(/^\#/);
    my ($chr, $start, $end, $id, $expr, $strand) = split(/\s+/, $_);
    if($expr == 0){$expr = 1;}
    if(!exists($unique{$id})){$readCount += $expr;}    
    $unique{$id}++;
    $length{($end-$start+1)}+=$expr;
    $entryCount++;
  }
  close(FILE);

  # WRITE LENGTH DISTRIBUTION
  open(LEN, ">$lengthDistrFile") || die "cannot open $lengthDistrFile\n";
  foreach my $length (sort {$a<=>$b} keys %length){
    print LEN "$length\t$length{$length}\n";
  }
  close(LEN);
  %length = ();

  # WRITE MULTIPLE MAPPINGS DISTRIBUTION
  my %multi = (); my $uniqueMappingReads = 0;
  foreach my $id (keys %unique){
    $multi{$unique{$id}}++;
  }
  open(MULTI, ">$multiMapFile") || die "cannot open $multiMapFile\n";
  foreach my $m (sort {$a<=>$b} keys %multi){
    if($m == 1){$uniqueMappingReads = $multi{$m};}
    print MULTI "$m\t$multi{$m}\n"
  }
  close(MULTI);
  %multi = ();

  # WRITE INFO
  open(OUT, ">$summaryFile") || die "cannot open $summaryFile\n";
  my $readCount_ts = $readCount;
  $tagCount = keys(%unique);
  my $tagCount_ts = $tagCount;
  my $entries_ts = $entryCount;
  my $uniqueMappingReads_ts = $uniqueMappingReads;
  1 while $readCount_ts =~ s/^(-?\d+)(\d{3})/$1,$2/;
  1 while $tagCount_ts =~ s/^(-?\d+)(\d{3})/$1,$2/;
  1 while $entries_ts =~ s/^(-?\d+)(\d{3})/$1,$2/;
  1 while $uniqueMappingReads_ts =~ s/^(-?\d+)(\d{3})/$1,$2/;

  print OUT "reads:$readCount:$readCount_ts\ntags:$tagCount:$tagCount_ts\nentries:$entryCount:$entries_ts\nuniqueMappingTags:$uniqueMappingReads:$uniqueMappingReads_ts\n";
  close(OUT);
    
  return 1;
}

sub writeNormalizedBed{
    my ($file, $outFile) = @_;
    my $id = 0; my $currentChrom = "NA";
    my %tags = ();
    
    
    my %chroms = ();
    open(FILE, "<$file") || die "cannot open $file\n";
    while(<FILE>){
        chomp;
        next if(/^\#/);
        my ($chr, $start, $end, $id, $expr, $strand) = split(/\s+/, $_);
        
        # if chromosome is written incorrect -> change it
        $chr =~ s/chrom/chr/;
        $chr =~ s/chromosome/chr/;
        if($chr =~ /^\d+$/){$chr = "chr".$chr;}
        if($chr =~ /^X$/){$chr = "chr".$chr;}
        if($chr =~ /^Y$/){$chr = "chr".$chr;}
        $chroms{$chr} = 1;     
    }
    close(FILE);
    
    foreach my $chrom (keys %chroms){
        my @lines = `cat $file | grep $chrom`;
        #        open(FILE, "<$file") || die "cannot open $file\n";
        #        while(<FILE>){
        #            chomp;
        foreach (@lines){
            next if(/^\#/);
            my ($chr, $start, $end, $id, $expr, $strand) = split(/\s+/, $_);
            if($expr == 0){$expr = 1;}
            
            # if chromosome is written incorrect -> change it
            $chr =~ s/chrom/chr/;
            $chr =~ s/chromosome/chr/;
            if($chr =~ /^\d+$/){$chr = "chr".$chr;}
            if($chr =~ /^X$/){$chr = "chr".$chr;}
            if($chr =~ /^Y$/){$chr = "chr".$chr;}
            if($chr =~ /^M$/){$chr = "chr".$chr;}
            
            next if($chr ne $chrom);
            
            $tags{$chr}{$strand}{$start}{$end}{expr} += ($expr / $unique{$id});
            $tags{$chr}{$strand}{$start}{$end}{readCnt} += $expr;
        }
        $id = printResult($outFile, $id, %tags);
        %tags = ();
        close(FILE);
    }
}

sub printResult{
    my ($outFile, $id, %tags) = @_;      
    open(OUT, ">>$outFile") || die "cannot open $outFile\n";
    foreach my $chrom (sort {$a cmp $b} keys %tags){
        foreach my $strand (keys %{$tags{$chrom}}){
            foreach my $start (sort {$a <=> $b} keys %{$tags{$chrom}{$strand}}){
                foreach my $end (sort {$a <=> $b} keys %{$tags{$chrom}{$strand}{$start}}){
                    $id++;
                    print OUT "$chrom\t$start\t$end\tdario\_$id\t$tags{$chrom}{$strand}{$start}{$end}{expr}\t$strand\t$tags{$chrom}{$strand}{$start}{$end}{readCnt}\n";
                }
            }
        }
    }
    close(OUT);
    return $id;
}


