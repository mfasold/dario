#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Getopt::Std;
use List::Util;
use Cwd;
use IO::File;
use POSIX qw(tmpnam);

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $inputFile $outputFile $format $zip $readFilter);
$readFilter = 1;

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions (
"i=s"       => \$inputFile,
"o=s"       => \$outputFile,
"f=s"       => \$format,
"z"         => \$zip,
"r=s"         => \$readFilter,
"help"      => \$help,
"h"         => \$help);
usage() if ($help || !$inputFile || !$outputFile || !($format == 1 || $format == 2 || $format == 3));


# -----------------------------------------------------------------------------
# MAIN

printHeader();

if($format == 1 || $format == 3){
    sam2bed($inputFile, $outputFile);
}
if($format == 2){
    soap2bed($inputFile, $outputFile);
}

if(!$zip){
    my $status = system("gzip $outputFile");
    if($status == 0){print "generated output file: $outputFile\.gz\n\n";}
    else{print "\nmap2bed.pl:\ngzip not installed -> Please pack you output file manually, to avoid huge upload files!\n\ngenerated output file: $outputFile\n\n";}
}


# -----------------------------------------------------------------------------
# FUNCTIONS

sub usage {
    print STDERR "\nusage: map2bed.pl -i <file> -f <int> -o <file>\n";
    print STDERR "bring mapping output to bed format\n";
    print STDERR "\n";
    print STDERR "[INPUT]\n";
    print STDERR " -i <file>    mapped reads file\n";
    print STDERR " -f <int>     format:\n";
    print STDERR "               1: SAM\n";
    print STDERR "               2: SOAP\n";
    print STDERR "               3: BAM\n";
    print STDERR " -o <file>    output file\n";
    print STDERR " -r <int>     read filter (minimum number of reads a\n";
    print STDERR "              tag must have to be written) (default: 1)\n";
    print STDERR " -h           this (usefull) help message\n";
    print STDERR "[VERSION]\n";
    print STDERR " 07-02-2013\n";
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

sub printHeader{
    print "# map2bed.pl started " . prettyTime() . "\n";
    print "# input file: $inputFile\n";
    print "# file format: ";
    if($format == 1){print "SAM\n";}
    if($format == 2){print "SOAP\n";}
    if($format == 3){print "BAM\n";}
    print "#\n";
}

sub checkFormat{
    my ($entry, $format) = @_;
    if($format eq "SOAP"){
        my @l = split(/\t+/, $entry);
        if($l[8] !~ /\d+/ || ($l[6] ne "+" && $l[6] ne "-")){die "file is not in SOAP format -> $entry\n";}
    }
    if($format eq "SAM"){
        my @l = split(/\t+/, $entry);
        if($l[3] !~ /\d+/){die "file is not in SAM format -> $entry\n";}
        if($l[1]&64 || $l[1]&128){die "It seems your file consists of paired end reads.\nUnfortunately, we do not support these experiments.\n";}
    }
}


sub sam2bed{
    my ($inputFile, $outputFile) = @_;
    my %tags = (); my $c = 0;
    if($format == 1){
        open(FILE, "<$inputFile") || die "cannot open $inputFile\n";
    }else{
        open(FILE, "samtools view $inputFile | ") || die "cannot open $inputFile\n";
    }
    while(<FILE>){	
        chomp;
        my @entry = split(/\t/,$_);
        next if ($_ =~ /^[@]/ || $_ =~ /^\s*$/ || $entry[2] eq "*");        
        checkFormat($_, "SAM");
        my ($id, $chr, $strand, $seq, $start, $end, $cigar) = sam($_);
        
        # name tag
        if(!exists($tags{$seq}{tag})){
            $c++;
            $tags{$seq}{tag} = "tag_$c";            
        }
        # count reads
        if(!exists($tags{$seq}{reads}{$id})){
            $tags{$seq}{reads}{$id} = 1;
        }
        # count loci
        if($cigar !~ /N/){
            $tags{$seq}{loci}{"$chr|$start|$end|$strand"} = 1;
        }else{
            my @parts = getSpliceParts($cigar, $start);
            foreach my $part (@parts){
                ($start, $end) = split(/\|/, $part);
                $tags{$seq}{loci}{"$chr|$start|$end|$strand"} = 1;
            }
        }
        
        foreach my $col (@entry){
            if($col =~ /XA:Z:/ && $col !~ /XA:Z:Q/){
                $col =~ s/XA:Z://;
                my @multi = split(/\;/,$col);
                foreach my $mm (@multi){
                    ($chr, $start, $cigar) = split(/\,/, $mm);
                    if($start =~ /\+/){$strand = "+";}else{$strand = "-";}
                    $start =~ s/[\+\-]//;
                    my $length = 1; my $tmp = $cigar; $tmp =~ s/(\d+)[MD]/$length+=$1/eg;
                    my $end = $start + $length;
                    if($cigar !~ /N/){
                        $tags{$seq}{loci}{"$chr|$start|$end|$strand"} = 1;
                    }
                    else{
                        my @parts = getSpliceParts($cigar, $start);
                        foreach my $part (@parts){
                            ($start, $end) = split(/\|/, $part);
                            $tags{$seq}{loci}{"$chr|$start|$end|$strand"} = 1;
                        }
                    }
                }
            }
        }
    }
    close(FILE);
    
    open(OUT, ">$outputFile") || die "cannot open $outputFile\n";
    foreach my $seq (keys %tags){
        foreach my $locus (keys %{$tags{$seq}{loci}}){
            my ($chr, $start, $end, $strand) = split(/\|/,$locus);
            if(keys(%{$tags{$seq}{reads}}) >= $readFilter){
                print OUT "$chr\t$start\t$end\t$tags{$seq}{tag}\t".(keys (%{$tags{$seq}{reads}}))."\t$strand\t".(keys (%{$tags{$seq}{loci}}))."\n";
            }
        }
    }
    close(OUT);
}

sub sam{
    my ($line) = @_;
    my @entry = split(/\t/,$line);
    my $id = $entry[0];
    my $chr = $entry[2];
    my $strand = $entry[1]?"-":"+";
    my $seq = $entry[9]; if($strand eq "-"){$seq = reverseComplement($seq);}
    my $start = $entry[3];
    my $cigar = $entry[5];
    my $length = 1; my $tmp = $cigar; $tmp =~ s/(\d+)[MD]/$length+=$1/eg;
    my $end = $entry[3] + $length;
    return ($id, $chr, $strand, $seq, $start, $end, $cigar);
}

sub soap2bed{
    my ($inputFile, $outputFile) = @_;
    my %bins = ();
    
    # split sam file into x bins (x = nb of chroms)
    open(FILE, "<$inputFile") || die "cannot open $inputFile\n";
    while(<FILE>){	
        chomp;
        next if (/^\#/);
            checkFormat($_, "SOAP");
        my @entry = split(/\s+/, $_);
        my $chr = $entry[7]; if($chr !~ /chr/){$chr = "chr".$chr;}
        my $strand = $entry[6];
        my $start = $entry[8];
        my $end = $entry[8] + length($entry[1]);  
        if(!exists($bins{$chr})){    
            ($bins{$chr}{fileName}, $bins{$chr}{fileHandle}) = openBin();
        }
        my $file = $bins{$chr}{fileHandle};
        print $file "$chr|$start|$end|$entry[1]|$strand";
    }
    
    close(FILE); 
    
    # run through the bins and merge the reads to tags
    my $c = 0;
    open(OUT, ">$outputFile") || die "cannot open $outputFile\n";
    foreach my $chrom (keys %bins){
        my %hash = ();
        open(TMP, "<$bins{$chrom}{fileName}") || die "cannot open $chrom->$bins{$chrom}{fileName}\n";
        while(<TMP>){
            chomp;
            $hash{$_}++;
        }
        close(TMP);
        foreach my $entry (keys %hash){
            my @bed_a = split(/\|/, $entry);
            next if(!$bed_a[4]);
            $c++;
            $bed_a[3] = "tag_$c";
            $bed_a[5] = $bed_a[4];
            $bed_a[4] = $hash{$entry};
            print OUT join("\t", @bed_a)."\n";
        }
    }
    close(OUT);
    closeBins(%bins);
} 

sub getSpliceParts{
    my ($cigar, $start) = @_;
    my @parts = ();
    my $length = 0; my $value = ""; 
    for(my $i = 0; $i <= length($cigar); $i++){
        my $v = substr($cigar,$i,1);
        if($v=~/\d/){$value.=$v;}
        if($v=~/\D/){
            if($v eq "M" || $v eq "D"){$length+=$value;$value="";}
            if($v eq "N"){
                my $end += $start + $length + 1;
                push(@parts, "$start|$end");
                $start = $start + $length + $value;
                $length = 0;$value="";}}}
    my $end += $start + $length + 1;
    push(@parts, "$start|$end");
    return @parts;
}

sub reverseComplement{
    my ($seq) = @_;
    my $revcomp = reverse($seq);    
    $revcomp =~ tr/ACGTUacgtu/TGCAAtgcaa/;    
    return $revcomp;
}
