#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);

my %reruns;

# chunks that should rerun
my $rerunF = '/home/dsellis/data/IES/analysis/allvsall/2rerun.dat';

# results from cluster
my $blastoutD = '/home/dsellis/data/IES/analysis/allvsall/blastout/';

#read files that need to be rerun

open IN, $rerunF or die $!;
while(my $line = <IN>){
    $line =~ /^(\d+),.*/;
    $reruns{$1} = 1;
}
close IN;

opendir(DH, $blastoutD) or die $!;
my @nodes = readdir(DH);
close DH;
my %clusterF;
foreach my $node (@nodes){
    next if $node eq '.' or $node eq '..';
    my $nodePath = catfile($blastoutD, $node);
    if (-d $nodePath){
	print $nodePath,"\n";
	opendir(ND, $nodePath);
	my @files = grep {/\.dat$/} readdir(ND);
	print "@files\n";
	foreach my $file (@files){
	    $file =~ /blastout\.chunk\.(\d+)\.fa\.dat/;
	    my $chunk = $1;
	    if(defined($clusterF{$file})){
		print $file," duplicated \n";
	    }else{
		if(defined($reruns{$chunk})){
		    print $file," should rerun\n";
		}else{
		    $clusterF{$file} = $node
		}
	    }
	}
	close ND;
    }
}

# gather blast results from cluster, find missing ones and run locally
# check if there are duplicate runs
# if not move to common folder
# rerun locally
