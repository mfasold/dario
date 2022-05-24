#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $bedFile $outFolder $annotationFile %unique %tags %entries $srcDir);

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("i=s"       => \$bedFile,
            "u=s"       => \$annotationFile,
            "o=s"       => \$outFolder,
            "w=s"       => \$srcDir,
            "help"      => \$help,
            "h"         => \$help);
usage() if ($help || !$bedFile || !$outFolder);

open(ERROR, ">>$outFolder/run.log") || die "cannot open $outFolder/run.log\n";
print ERROR "checkUpload.pl: started (".prettyTime().")\n";



# -----------------------------------------------------------------------------
# MAIN

my ($bedFile, $ff) =  extract($bedFile);
print ERROR "checkUpload.pl: mapping file was $ff packed\n";

if($annotationFile){
  my ($annotationFile, $ff) =  extract($annotationFile);
  print ERROR "checkUpload.pl: annotation file was $ff packed\n";
}

($bedFile, $ff) =  checkFileFormat($bedFile);
print ERROR "checkUpload.pl: mapping file was in $ff format\n";

my $result = check($bedFile, "upload");
if($result == 0){
  print ERROR "checkUpload.pl: mapping file has wrong file format\n";
  die "mapping file has wrong file format";
}else{
  writeNormalizedBed($bedFile, "$bedFile\.tmp");
  system("gzip $bedFile");
  system("mv $bedFile\.tmp $bedFile");
  print ERROR "checkUpload.pl: mapping bed format correct\n";
}


if($annotationFile){
  my $result = check($annotationFile, "annotation");
  if($result == 0){
    print ERROR "checkUpload.pl: annotation file has wrong file format\n";
    die "user_annotation file has wrong file format";
  }else{
    print ERROR "checkUpload.pl: annotation bed format correct\n";
  }
}

print ERROR "checkUpload.pl: done (".prettyTime().")\n";
close(ERROR);




# -----------------------------------------------------------------------------
# FUNCTIONS

sub usage {
  print STDERR "\nusage: checkUpload.pl -i <file> -o <dir>\n";
  print STDERR "check bed file\n";
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

sub check{
  my ($file, $flag) = @_;
  my $line = 0;
  my %length = ();
  my $readCount = 0;
  my $tagCount = 0;
  my $entryCount = 0;

  # FIND MULTIPLE MAPPINGS
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    $line++;
    next if(/\#/);
    my ($chr, $start, $end, $id, $expr, $strand) = split(/\s+/, $_);
    # check for correct entry types
    if($start !~ /^\d+$/ || $end !~ /^\d+$/ || $expr !~ /^\d+$/ || ($strand ne "+" && $strand ne "-")){print ERROR "checkUpload.pl: error in line $line -> \"$_\"\n"; return 0;}
    if(!exists($unique{$id})){$readCount += $expr;}    
    $unique{$id}++;
    $length{($end-$start+1)}+=$expr;
    $entryCount++;
  }
  close(FILE);

  if($flag eq "annotation"){return 1;}

  # WRITE LENGTH DISTRIBUTION
  open(LEN, ">$outFolder/length.out") || die "cannot open $outFolder/length.out\n";
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
  open(MULTI, ">$outFolder/multipleMappings.out") || die "cannot open $outFolder/multipleMappings.out\n";
  foreach my $m (sort {$a<=>$b} keys %multi){
    if($m == 1){$uniqueMappingReads = $multi{$m};}
    print MULTI "$m\t$multi{$m}\n"
  }
  close(MULTI);
  %multi = ();

  # WRITE INFO
  open(OUT, ">$outFolder/upload.info") || die "cannot open $outFolder/upload.info\n";
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

sub extract{
    # TODO: generalize, do not use hard-coded "upload.bam" 
    my ($file) = @_;
    if($file =~ /\.tar\.gz$/){
        system("tar -C $outFolder -xvf $file");
        $file =~ s/\.tar\.gz//;
        # make sure the uncompressed files name will be upload.bed/upload.bam 
        system("find $outFolder -name '*.bed' -print0 | xargs -0 -I file mv file $outFolder/upload.bed");
        system("find $outFolder -name '*.bam' -print0 | xargs -0 -I file mv file $outFolder/upload.bam");
        return ($file, ".tar.gz");
    }
    elsif($file =~ /\.gz$/){
        system("gunzip $file"); # gzip automatically unzips in source folder and names extraxt upload.*
        $file =~ s/\.gz//;
        return ($file,".gz");
    }
    elsif($file =~ /\.zip/){
        system("unzip -d $outFolder $file");
        $file =~ s/\.zip//;
        # make sure the uncompressed files name will be upload.bed/upload.bam 
        system("find $outFolder -name '*.bed' -print0 | xargs -0 -I file mv file $outFolder/upload.bed");
        system("find $outFolder -name '*.bam' -print0 | xargs -0 -I file mv file $outFolder/upload.bam");
        return ($file,".zip");
    }
    else{
        return ($file,"not");
    }
}

sub checkFileFormat{
    my ($file) = @_;
    if($file =~ /\.bam$/){
        my $newFileName = $file;
        $newFileName =~ s/\.bam//;
        system("samtools view -h $file > $newFileName\.sam");
        system("$srcDir/analysis/map2bed.pl -i $newFileName\.sam -f 1 -o $newFileName -z");
        system("rm $newFileName\.sam");
        return($newFileName, "bam");
    }
    elsif($file =~ /\.bed$/){
        return($file, "bed");
    }
    else{
        die "Wrong file extension! Only .bed, or .bam!\n";
    }
}

