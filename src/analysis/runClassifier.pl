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

my ($testFile, $modelFile, $predictionsFile, $help, $outFolder);
GetOptions( "testSet=s"     => \$testFile,
	    "model=s"       => \$modelFile,
	    "predictions=s" => \$predictionsFile,
	    "o=s"           => \$outFolder,
	    "help"          => \$help);
if(!$testFile || !$modelFile || !$predictionsFile || !$outFolder || $help){usage(); exit;}

#open(ERROR, ">>$outFolder\/run.log") || die "cannot open $outFolder\/run.log\n";
print "runClassifier.pl: started (".prettyTime().")\n";

# read test set
my $minBlocks = 2;
my %testSet = readBlocks($testFile);

# calc features
my %features = getFeatures(%testSet);

# write arff file for weka
my %adr = writeFeatures("test.arff");

# run weka to create model
my @classification = `java -classpath \$CLASSPATH:$weka -Xmx256M weka.classifiers.trees.RandomForest -l $modelFile -T test.arff -p 0 -distribution`;
system("rm test.arff");

# parse results
open(OUT, ">$predictionsFile") || die "cannot open $predictionsFile\n";
my %statistics = parseResult(@classification);
close(OUT);

# write statistics file
foreach my $type (keys %statistics){
  if(!$statistics{$type}{false}){
    print "$type\t$statistics{$type}{true}/$statistics{$type}{true} (1)\n";
  }
  elsif(!$statistics{$type}{true}){
    print "$type\t0/$statistics{$type}{false}(0)\n";
  }
  else{
    print "$type\t$statistics{$type}{true}/".($statistics{$type}{true} + $statistics{$type}{false})." (".($statistics{$type}{true} / ($statistics{$type}{true} + $statistics{$type}{false})).")\n";
  }
}



print "runClassifier.pl: done (".prettyTime().")\n";
#close(ERROR);


#############
# functions #
#############

sub usage{
    print "usage: ./trainClassifier.pl --trainSet file <FILE> --model <FILE> --statistic <FILE>\n";
    print "\n";
    print "[INPUT]\n";
    print " -trainSet <FILE>\tfile with flagged clusters \n";
    print " -model <FILE>\tfile with the calculated model\n";
    print "[OUTPUT]\n";
    print " -predictions <FILE>\tfile with the new predictions\n";
    print " -help\tthis (helpfull) message\n";
    print "[VERSION]\n";
    print " 07-21-2010\n";
    print "[BUGS]\n";
    print " Please report bugs to david\@bioinf.uni-leipzig.de\n";
    print "\n";
}

sub printHeader{
  my ($file) = @_;
  open(STATISTICS, ">$predictionsFile") || die "cannot open $predictionsFile\n";
  print STATISTICS "# trainClassifier.pl output generated " . prettyTime() . "\n";
  print STATISTICS "# training set: $testFile\n";
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
    next if(/\#/);
    if(/\>/){
      my ($clusterID, $chrom, $start, $end, $strand, $absHeight, $relHeight, $blockCount) = split(/\s+/, $_);
      next if($blockCount < $minBlocks);
      my @check = split(/\s+/, $_);

      $hash{$clusterID}{chrom} = $chrom;
      $hash{$clusterID}{start} = $start;
      $hash{$clusterID}{end} = $end;
      $hash{$clusterID}{strand} = $strand;
      $hash{$clusterID}{blockCount} = $blockCount;
      $hash{$clusterID}{absHeight}  = $absHeight;
      $hash{$clusterID}{relHeight}  = $relHeight;
      
      if($check[14]){
	if($check[14] eq "miRNA" || $check[14] eq "tRNA" || $check[14] eq "snoRNA_HACA" || $check[14] eq "snoRNA_CD"){$hash{$clusterID}{ncRNAtype} = $check[14];}
	else{$hash{$clusterID}{ncRNAtype} = "ncRNA";}
      }
      elsif($check[8]){$hash{$clusterID}{ncRNAtype} = "ncRNA";}
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
  my %hash = (); my $c = 0;
  open(FILE, ">$file") || die "cannot open $file\n";

  print FILE "\@RELATION HTS_DATA\n";
  foreach my $feature (@myFeatures){
    print FILE "\@ATTRIBUTE $feature\tNUMERIC\n";
  }
  print FILE "\@ATTRIBUTE class\t\{miRNA, snoRNA_HACA, snoRNA_CD, tRNA\}\n";
  print FILE "\@DATA\n";
  foreach my $id (keys %features){
    $c++;
    foreach my $feature (@myFeatures){
      print FILE "$features{$id}{$feature},";
    }	    
    if($testSet{$id}{ncRNAtype} eq "tRNA" || $testSet{$id}{ncRNAtype} eq "snoRNA_HACA" || $testSet{$id}{ncRNAtype} eq "snoRNA_CD"  || $testSet{$id}{ncRNAtype} eq "miRNA"){
      print FILE "$testSet{$id}{ncRNAtype}\n";
    }
    else{
      print FILE "?\n";
    }
    $hash{$c} = $id;
  }
  close(FILE);
  return %hash;
}

sub parseResult{
  my (@classification) = @_;
  my %statistics = (); my $c = 0;

  foreach (@classification){
    chomp;
    $_ =~ s/\*//;
    $_ =~ s/\+//;
    my ($dum, $inst, $actual, $predicted, $dist) = split(/\s+/, $_);
    next if($_ =~ /^$/ || $_ =~ /^\s+$/);
    next if($inst !~ /\d+/);
    my ($microRNA, $snoRNA_HACA, $snoRNA_CD, $tRNA) = split(/\,/, $dist);
    
    if($actual =~ /miRNA/ && $predicted =~ /miRNA/){$statistics{miRNA}{true}++;}
    if($actual =~ /miRNA/ && ($predicted =~ /tRNA/ || $predicted =~ /snoRNA/)){$statistics{miRNA}{false}++;}

    if($actual =~ /tRNA/ && $predicted =~ /tRNA/){$statistics{tRNA}{true}++;}
    if($actual =~ /tRNA/ && ($predicted =~ /miRNA/ || $predicted =~ /snoRNA/)){$statistics{tRNA}{false}++;}

    if($actual =~ /snoRNA/ && $predicted =~ /snoRNA/){$statistics{snoRNA}{true}++;}
    if($actual =~ /snoRNA/ && ($predicted =~ /tRNA/ || $predicted =~ /miRNA/)){$statistics{snoRNA}{false}++;}    

    if($testSet{$adr{$inst}}{ncRNAtype} eq "NA"){
      my $score = 0; my $type = "";
      if($predicted =~ /miRNA/){$score = $microRNA; $type = "miRNA";}
      if($predicted =~ /tRNA/){$score = $tRNA; $type = "tRNA";}
      if($predicted =~ /snoRNA_C/){$score = $snoRNA_CD; $type = "snoRNA_CD";}
      if($predicted =~ /snoRNA_H/){$score = $snoRNA_HACA; $type = "snoRNA_HACA";}
      $c++;
      my $string = "$testSet{$adr{$inst}}{chrom}\t$testSet{$adr{$inst}}{start}\t$testSet{$adr{$inst}}{end}\t$type\_$c\t$score\t\t$testSet{$adr{$inst}}{strand}\t$type";
      my $folderID = $outFolder; $folderID =~ s/\/scratch\/dario\/computations\///; $folderID =~ s/\///;
      $string =~ s/\"//g;
      print OUT $string."\n";

    }
  }
  return %statistics;
}
