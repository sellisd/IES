#!/usr/bin/perl
use warnings;
use strict;
#parse ancestral state reconstructions from rev bayes and prepare it for easy processing
my $asrPath = '/home/dsellis/data/IES_data/msas/asr/';
my @runs = qw/run4 run5/;
open GF, '/home/dsellis/data/IES_data/msas/asr/geneFamilies.dat' or die;
my @clusters = <GF>;
close GF;
map(chomp, @clusters);

my $burnIn = 1000;
print "cluster\tnode\tiesColumn\tpresence\n";
foreach my $cluster (@clusters){
    my @nodes;
    my $iesColumns = 0;
    my @sums;
    my $its = 0;
    foreach my $run (@runs){
	my $file = $asrPath.$run.'/ancStates'.$cluster.'.log';
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
	    print join("\t", ($cluster, $node, $col, $sums[$counter]/$its)), "\n";
	    $counter++;
	}
    }
    close IN;
}
