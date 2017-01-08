#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Spec::Functions qw(catfile);
use lib'.';
use functions;

my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};
my $help;
my $asrPath;
my $outputF;
my $gf;
my $burnIn = 1000;
my $usage = <<HERE;

Parse ancestral state reconstructions from revBayes output and prepare it in a tidy format
usage parseRevBayesAsr.pl [OPTIONS]
where OPTIONS can be:
  -asr:    path for one or more revbayes runs (each in a runX file, where X = 1,2,...)
  -output: output file
  -burnin: burn in period
  -gf:     file with list of gene families expected
  -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?'   => \$help,
			      'output=s' => \$outputF,
			      'burnin=s' => \$burnIn,
			      'asr=s'    => \$asrPath,
			      'gf=s'     => \$gf
		   ));

die $usage if $help;


my @runs = qw/run1 run2/;
open GF, $gf or die;
my @clusters = <GF>;
close GF;
map(chomp, @clusters);

open OUT, '>', $outputF or die $!;
print OUT "cluster\tnode\tiesColumn\tpresence\n";
foreach my $cluster (@clusters){
    my @nodes;
    my $iesColumns = 0;
    my @sums;
    my $its = 0;
    foreach my $run (@runs){
	my $file = catfile($asrPath, $run, 'ancStates'.$cluster.'.log');
	if (! -e $file){
	    print "skipping $file\n";
	    next;
	}
	open IN, $file or die "$! $file";
	my $lineCounter = 0;
	while(my $line = <IN>){
	    chomp $line;
	    if($lineCounter == 0 and $run eq $runs[0]){ # build header from first file processed
		#read header
		my @header = split " ", $line;
		shift @header; # iteration
		foreach my $node (@header){
		    $node =~ s/end_//;
		    push @nodes, $node;
		}
	    }elsif($lineCounter > 0){ # exclude header to process all files
		my @ar = split " ", $line;
		next unless $ar[0] > $burnIn;
		if($lineCounter == 1 and $run eq $runs[0]){ # find number of ies columns (first file only)
		    my @gauge = split ',', $ar[1];
		    $iesColumns = $#gauge + 1;
		    # initialize array of colsums only for the first file processed
		    if($run eq $runs[0]){
			for(my $i = 0; $i < $iesColumns * ($#nodes + 1); $i++){
			    $sums[$i] = 0;
			}
		    }
		}
		shift @ar; # remove Iterations
		my @all;
		foreach my $str (@ar){
		    my @pa = split ",", $str;
		    push @all, @pa;
		}
		for(my $i = 0; $i <= $#all; $i++){
		    $sums[$i] += $all[$i];
		}
		$its++;
	    }
	    $lineCounter++;
	}
    }
    my $counter = 0;
    foreach my $node (@nodes){
	for (my $col = 0; $col < $iesColumns; $col++){
	    print OUT join("\t", ($cluster, $node, $col, $sums[$counter]/$its)), "\n";
	    $counter++;
	}
    }
    close IN;
}
close OUT;
