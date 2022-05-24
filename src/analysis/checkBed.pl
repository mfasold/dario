#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $bedFile);

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("i=s"       => \$bedFile,
            "help"      => \$help,
            "h"         => \$help);
usage() if ($help || !$bedFile);


# -----------------------------------------------------------------------------
# MAIN
print "checkBed.pl: check uploaded bed file (".prettyTime().")\n";

checkBed($bedFile);

print "checkBed.pl: .bed format correct (".prettyTime().")\n";


# -----------------------------------------------------------------------------
# FUNCTIONS

sub usage {
  print STDERR "\nusage: checkBed.pl -i <file>\n";
  print STDERR "check for correct bed format\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -i <file>    bed file\n";
  print STDERR " -h <file>    this (usefull) help message\n";
  print STDERR "[VERSION]\n";
  print STDERR " 06-07-2010\n";
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

sub checkBed{
  my ($file) = @_;
  my %unique = ();
  my %length = ();
  my $line = 0;
  
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    $line++;
    next if(/^\#/); # run over header
        if(/^track/ && $line == 1){system("sed -i 1d $bedFile"); print "Galaxy header deleted\n"; next;}
    my @array = split(/\s+/, $_); 
      #print $#array."\n$array[0]\t$array[5]\n";
      if($#array < 5){
       printError("bed", $line, $_, "No Header and not enough columns."); die "checkBed.pl:\nerror in line $line -> \"$_\"\nNot enough columns.\n";    
      } 
    my ($chr, $start, $end, $id, $expr, $strand) = split(/\s+/, $_); # split in correct bed format
      
# read to long for short RNA-seq protocol
      if((($end - $start) + 1) > 55){
          print "Your dataset comprises mapping loci longer than 55nt.\n";
          print "For now, we only support reads produced using a short RNA-seq protocol.\n"; 
          die "checkBed.pl: not a short RNA-seq protocol.\n";
      }
      
# correct wrong chromosome names
      $chr =  lc($chr);
      $chr =~ s/chromosome/chr/;
      $chr =~ s/chrom/chr/;
      if($chr =~ /^\d+$/){$chr = "chr".$chr;}
   
# chromosome is not in correct format
#      if($chr  !~ /^chr\S+$/ ){
#          printError("bed", $line, $_, "Chromosome is not in correct format."); 
#          die "checkBed.pl:\nerror in line $line -> \"$_\"\nChromosome is not in correct format.\n";
#      }
# start is not numeric
      if($start !~ /^\d+$/){
          printError("bed", $line, $_, "Chromosome is not in correct format."); 
          die "checkBed.pl:\nerror in line $line -> \"$_\"\nStart position is not numeric.\n";
      } 
# end is not numeric
      if($end !~ /^\d+$/){
          printError("bed", $line, $_, "End position is not numeric."); 
          die "checkBed.pl:\nerror in line $line -> \"$_\"\nEnd position is not numeric.\n";
      } 
# expression is not numeric
      if($expr !~ /^\d+$/){
          printError("bed", $line, $_, "Expression value is not numeric."); 
          die "checkBed.pl:\nerror in line $line -> \"$_\"\nExpression value is not numeric.\n";
      } 
# expression cannot be 0
#         if($expr == 0){
#          printError("bed", $line, $_, "Expression cannot be zero."); 
#          die "checkBed.pl:\nerror in line $line -> \"$_\"\nExpression cannot be zero.\n";
#      } 
# strand is not + or -
      if($strand ne "+" && $strand ne "-"){
          printError("bed", $line, $_, "Strand is not + or -."); 
          die "checkBed.pl:\nerror in line $line -> \"$_\"\nStrand is not + or -.\n";
      }

      $unique{$id}++;
      my $l = ($end-$start)+1;$length{$l}++;
  }
# IDs are not unique
  if(keys(%unique) < 2){
    print "IDs seem not to be unique for different reads/tags\n"; 
    die "checkBed.pl: IDs seem not to be unique for different reads/tags\n";}
# IDs are not unique
  if(keys(%length) < 2){
    print "All mapping loci have the same length.\nPlease check your dataset and/or mapping tool.\n"; 
    die "checkBed.pl: All mapping loci have the same length.\n";}

  close(FILE);
}

sub printError{
    my ($type, $nb, $line, $message) = @_;
    
    if($type eq "bed"){
        print "Wrong file Format.\n\n";
        print "Line $nb of your file generates an error.\n";
        print "line $nb: $line\n\n";
        print "error message:\n$message\n\n";
        print "Please use the .bed format as follows:\n";
        print "#chr    start   end     id              nb of sequenced reads      strand\n";
        print "chr1    20229   20366   GSM34290_325    50                         +\n";
        print "chr1    20230   20369   GSM34290_328    13                         +\n";
        print "\n";
        print "This might help you to fix the problem:\n";
        print "- flag header with a \"\#\" symbol at the start of the line\n";
        print "- the columns in the bed-file have to be separated with tabs, no\n";
        print "  whitespaces are allowed\n";
        print "- all six columns have to be filled with values\n";
        print "- due to the overlapping we perform, the chromosome names have to\n";
        print "  be in the UCSC-style (chr1, chr2, chr3, ...)\n";
        print "- use different (unique) IDs for different reads/tags\n";
        print "- expression value must be greater than zero\n";
        print "- the strand has to be + or -\n";
    }
}










