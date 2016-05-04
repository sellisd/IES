#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# calculate coverage over regions defined in a bed file given a genome coverage file

# read bed file and make counters for each region
# read coverage file and fill in counters

my $bedF = '/home/dsellis/data/IES/analysis/bed/pca.ies.slop.sortmerg.be';
my $covF = '/home/dsellis/data/IES/genomicData/caudatum/genome/pca.cov';
my $outF = '/home/dsellis/data/IES/genomicData/caudatum/genome/pca.ies.cov';

print "reading bed file...";
open BED, $bedF or die $!;
my %regions;
my $curBegin = 0;
my $curChrom;
while(my $line = <BED>){
    chomp $line;
    (my $chrom, my $begin, my $end, my $name) = split "\t", $line;
    if(!defined($curChrom)){
	$curChrom = $chrom; # first instance
    }
    if($curChrom ne $chrom){
	$curChrom = $chrom;
    }else{
	if($curBegin > $begin){
	    die "bed file not sorted: $curBegin before $begin";
	}
    }
    $curBegin = $begin;
    if(defined($regions{$chrom})){
	push @{$regions{$chrom}}, {'begin' => $begin,
				   'end'   => $end,
				   'name'  => $name,
				   'cov'   => 0};
    }else{
	$regions{$chrom} = [{'begin' => $begin,
			     'end'   => $end,
			     'name'  => $name, 
			     'cov'   => 0}];
    }
}
close BED;
print "done\n";

print "reading coverage file\n";
my %chromCounter;
my $total = keys %regions;
open COV, $covF or die $!;
while(my $line = <COV>){
    chomp $line;
    (my $chr, my $pos, my $cov) = split "\t", $line;
    if(!defined($chromCounter{$chr})){
	$chromCounter{$chr} = 1;
	my $read = scalar keys %chromCounter;
	print $read, '/', $total,"\r";
    }
    for(my $i = 0; $i <= $#{$regions{$chr}}; $i++){
	# if current position is within region add coverage
	my $end = $regions{$chr}[$i]{'end'};
	my $begin = $regions{$chr}[$i]{'begin'};
	if($pos <= $begin){
	    last;
	}elsif($pos > $end){
	    next;
	}else{
	    $regions{$chr}[$i]{'cov'} += $cov;
	    last; # save time as we expect sorted
	}
	# if not continue
    }
}
close COV;
print "done\n";

open OUT, '>', $outF or die $!;
foreach my $chr (sort keys %regions){
    foreach my $be (@{$regions{$chr}}){
	my $covpernt = $be->{'cov'} / ($be->{'end'} - $be->{'begin'}); # per nt coverage
	print OUT join("\t", ($chr, $be->{'begin'}, $be->{'end'}, $covpernt, $be->{'name'})), "\n";
    }	
}
close OUT;
