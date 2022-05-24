#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;  
use IO::File;
use POSIX qw(tmpnam);

my $weka = "/scratch/dario/weka/weka.jar";

my @myFeatures = ("meanDist", 
		  "meanBlockLength", 
		  "maxDist", 
		  "blockCount", 
		  "minBlockLength", 
		  "maxHeight", 
		  "clusterLength", 
		  "maxBlockLength", 
		  "blockOverlapRange", 
		  "minHeight", 
		  "minDist", 
		  "blockOverlapHeight");

my ($trainFile, $modelFile, $statisticsFile, $help, $outFolder);
GetOptions( "trainSet=s"  => \$trainFile,
	    "model=s"     => \$modelFile,
	    "statistic=s" => \$statisticsFile,
	    "o=s"         => \$outFolder,
	    "help"        => \$help);
if(!$trainFile || !$modelFile || !$statisticsFile || !$outFolder || $help){usage(); exit;}

my $trainingFile = $outFolder."/train.arff";

# read training set
my $minBlocks = 2;
my %trainSet = readBlocks($trainFile);

# calc features
my %features = getFeatures(%trainSet);

# write arff file for weka
writeFeatures($trainingFile);

# run weka to create model
my @training = `java -classpath \$CLASSPATH:$weka -Xmx256M weka.classifiers.meta.AttributeSelectedClassifier -t $trainingFile -d $modelFile -E \"weka.attributeSelection.CfsSubsetEval\" -S \"weka.attributeSelection.BestFirst -D 1 -N 5\" -W weka.classifiers.trees.RandomForest -- -I 100 -K 0 -S 1`;
#system("rm train.arff");

# write statistics file
printHeader($statisticsFile);
open(STATISTICS, ">>$statisticsFile") || die "cannot open $statisticsFile\n";
foreach (@training){
  print "$_";
  print STATISTICS "$_";
}
close(STATISTICS);



#############
# functions #
#############

sub usage{
    print "usage: ./trainClassifier.pl --trainSet file <FILE> --model <FILE> --statistic <FILE>\n";
    print "\n";
    print "[INPUT]\n";
    print " -trainSet <FILE>\tfile with flagged clusters \n";
    print "[OUTPUT]\n";
    print " -model <FILE>\tfile with the calculated model\n";
    print " -statistic <FILE>\tfile with the model statistics\n";
    print " -help\tthis (helpfull) message\n";
    print "[VERSION]\n";
    print " 07-21-2010\n";
    print "[BUGS]\n";
    print " Please report bugs to david\@bioinf.uni-leipzig.de\n";
    print "\n";
}

sub printHeader{
  my ($file) = @_;
  open(STATISTICS, ">$statisticsFile") || die "cannot open $statisticsFile\n";
  print STATISTICS "# trainClassifier.pl output generated " . prettyTime() . "\n";
  print STATISTICS "# training set: $trainFile\n";
  print STATISTICS "#\n";
  close(STATISTICS);
}

sub prettyTime{
  my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
  my ($second, $minute, $hour, $dayOfMonth, $month, 
      $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  my $year = 1900 + $yearOffset;
  return "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
}

sub readBlocks{
  my ($file) = @_;
  my %hash = ();
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    if(/\>/){
      my ($clusterID, $chrom, $start, $end, $strand, $absHeight, $relHeight, $blockCount) = split(/\s+/, $_);
      next if($blockCount < $minBlocks);
      my @check = split(/\s+/, $_);
      
      $hash{$clusterID}{start} = $start;
      $hash{$clusterID}{end} = $end;
      $hash{$clusterID}{strand} = $strand;
      $hash{$clusterID}{blockCount} = $blockCount;
      $hash{$clusterID}{absHeight}  = $absHeight;
      $hash{$clusterID}{relHeight}  = $relHeight;
      
      if($check[9]){
	if($check[9] eq "miRNA" || $check[9] eq "tRNA" || $check[9] eq "snoRNA_HACA" || $check[9] eq "snoRNA_CD"){$hash{$clusterID}{ncRNAtype} = $check[9];}
	else{$hash{$clusterID}{ncRNAtype} = "NA";}
      }
      else{$hash{$clusterID}{ncRNAtype} = "NA";}
      
      for(my $i = 0; $i < $blockCount; $i++){
	my $block = <FILE>; chomp($block);
	my ($blockID, $blockChrom, $blockStart, $blockEnd, $blockStrand, $blockHeight) = split(/\s+/, $block);
	push(@{$hash{$clusterID}{normalizedBlockHeights}}, ($blockHeight / $absHeight));
	$hash{$clusterID}{blocks}{$blockID}{chrom}  = $blockChrom;
	$hash{$clusterID}{blocks}{$blockID}{start}  = $blockStart;
	$hash{$clusterID}{blocks}{$blockID}{end}    = $blockEnd;
	$hash{$clusterID}{blocks}{$blockID}{strand} = $blockStrand;
	$hash{$clusterID}{blocks}{$blockID}{height} = $blockHeight;
      }
    }
  }
  close(FILE);
  return %hash;
}

sub getFeatures{
    my (%cluster) = @_;
    my %features = ();
    foreach my $id (sort keys %cluster){
      my $dum;
      ($features{$id}{maxHeight}, $features{$id}{minHeight}, $dum) = parseArray(@{$cluster{$id}{normalizedBlockHeights}});
      ($features{$id}{maxDist}, $features{$id}{minDist}, $features{$id}{meanDist}) = getDist(%{$cluster{$id}{blocks}});
      $features{$id}{blockCount} = $cluster{$id}{blockCount};
#      ($features{$id}{maxBlockShift}, $features{$id}{minBlockShift}, $features{$id}{meanBlockShift}) = getBlockShift(%{$cluster{$id}{reads}});
      $features{$id}{clusterLength} = ($cluster{$id}{end} - $cluster{$id}{start});
      ($features{$id}{maxBlockLength}, $features{$id}{minBlockLength}, $features{$id}{meanBlockLength}) = getBlockLength(%{$cluster{$id}{blocks}});
      ($features{$id}{blockOverlapRange}, $features{$id}{blockOverlapHeight}) = getBlockOverlap($cluster{$id}{start}, $cluster{$id}{end}, %{$cluster{$id}{blocks}});
    }    
    return %features;
}

sub parseArray{
    my (@array) = @_;
    my ($max, $min, $mean);
    my $c = 0; my $sum = 0;
    $max = $array[0];
    $min = $array[0];
    foreach my $i (@array){
	$c++;
	$sum += $i;
        if ($i > $max){
            $max = $i;
        }
        elsif ($i < $min){
            $min = $i;
        }
    }
    $mean = $sum / $c;
    return($max, $min, $mean);
}

sub getDist{
    my (%blocks) = @_;
    my ($max, $min, $mean);
    my $last = -1; my @dists = ();
    foreach my $block (sort {$blocks{$a}{start} <=> $blocks{$b}{start} || $blocks{$a}{end} <=> $blocks{$b}{end}} keys %blocks){
	if($last == -1){
	    $last = $blocks{$block}{end};
	}
	else{
	    my $dist = $blocks{$block}{start} - $last;
	    push(@dists, $dist);
	    $last = $blocks{$block}{end};
	}
    }
    if($#dists == -1){return (0,0,0);}
    ($max, $min, $mean) = parseArray(@dists);
    return ($max, $min, $mean);
}

sub getBlockShift{
    my (%reads) = @_;
    my ($max, $min, $mean);
    my %hash = (); my @shifts = ();
    foreach my $read (keys %reads){
	$hash{$reads{$read}{block}}{height} += $reads{$read}{height};
	for(my $i = $reads{$read}{start}; $i <= $reads{$read}{end}; $i++){
	    $hash{$reads{$read}{block}}{$i}{height} += $reads{$read}{height};
	}
    }
    foreach my $block (keys %hash){
	my $length = 0; my $value = 0;
	foreach my $pos (keys %{$hash{$block}}){
	    next if($pos !~ /\d+/);
	    $length++;
	    $value += ($hash{$block}{$pos}{height} / $hash{$block}{height});
	}
	push(@shifts, ($value / $length));
    }
    ($max, $min, $mean) = parseArray(@shifts);
    return ($max, $min, $mean);
}

sub getBlockLength{
    my (%blocks) = @_;
    my ($max, $min, $mean);
    my @lengths = ();
    foreach my $block (keys %blocks){
	push(@lengths, ($blocks{$block}{end} - $blocks{$block}{start}));
    }
    ($max, $min, $mean) = parseArray(@lengths);
    return ($max, $min, $mean);

}

sub getBlockOverlap{
    my ($start, $end, %blocks) = @_;
    my ($max, $min, $mean);
    my @lengths = ();
    my $overlap = 0;
    my $height = 0;
    my $sum = 0;
    for(my $i = $start; $i <= $end; $i++){
	my $thisHeight = 0;
	foreach my $block (keys %blocks){
	    if($blocks{$block}{start} <= $i && $blocks{$block}{end} >= $i){
		$thisHeight++;
		if($thisHeight == 2){$overlap++;}
		if($thisHeight == 1){$sum++;}
	    }
	}
	if($thisHeight > $height){$height = $thisHeight;}
    }
    return (($overlap / $sum), $height);
}

sub writeFeatures{
  my($file) = @_;

  open(FILE, ">$file") || die "cannot open $file\n";

  print FILE "\@RELATION HTS_DATA\n";
  foreach my $feature (@myFeatures){
    print FILE "\@ATTRIBUTE $feature\tNUMERIC\n";
  }
  print FILE "\@ATTRIBUTE class\t\{miRNA, snoRNA_HACA, snoRNA_CD, tRNA\}\n";
  print FILE "\@DATA\n";
  foreach my $id (keys %features){
    if($trainSet{$id}{ncRNAtype} eq "tRNA" || $trainSet{$id}{ncRNAtype} eq "snoRNA_HACA" || $trainSet{$id}{ncRNAtype} eq "snoRNA_CD"  || $trainSet{$id}{ncRNAtype} eq "miRNA"){
      foreach my $feature (@myFeatures){
	print FILE "$features{$id}{$feature},";
      }	    
      print FILE "$trainSet{$id}{ncRNAtype}\n";
    }
  }
  close(FILE);
}
