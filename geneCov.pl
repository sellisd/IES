#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $help;
my $coverageF;
my $genesF;
my $out;
my $usage = <<HERE;

Calculate the average per nucleotide gene coverage over genomic regions. It reads the whole genome in memory, watch out for large genomes!

usage geneCov.pl [OPTIONS] INPUTFILE(S)
where OPTIONS can be:
  -cov:    coverage file
  -bed:    bed file with genomic regions (e.g genes)
  -out:    output file
  -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'cov=s'  => \$coverageF,
			      'bed=s'  => \$genesF,
			      'out=s'  => \$out
		   ));
die $usage if $help;

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
open OUT, '>', $out or die $!;
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
    print OUT $gene,"\t",$sum/($end - $begin),"\n";
}
close GN;
close OUT;
