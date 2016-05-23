#!/usr/bin/perl
use warnings;
use strict;

# calculate the average per base read coverage over each gene in *P. caudatum*

# reading whole genome in memory!

my $coverageF = '/home/dsellis/data/IES/genomicData/caudatum/genome/pca.cov';
my $genesF = '/home/dsellis/data/IES/analysis/bed/pca.gene.be';
my %h;

# read genomic coverage and make scaffold-sized arrays
open CV, $coverageF or die $!;
while(my $line = <CV>){
    chomp $line;
    (my $scaffold, my $location, my $coverage) = split " ", $line;
    ${$h{$scaffold}}[$location] = $coverage;
}
close CV;

# read genes, average slice from scaffold sized array

open GN, $genesF or die $!;
while(my $line = <GN>){
    chomp $line;
    (my $scaffold, my $begin, my $end, my $gene) = split " ", $line;
    my $sum = 0;
    for(my $i = ($begin+1); $i <= $end; $i++){
	if(defined($h{$scaffold}[$i])){
	    $sum += $h{$scaffold}[$i];
	} #if not defined sum is zero
    }
    print $gene,"\t",$sum/($end - $begin),"\n";
}
close GN;
