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

open COV, $covF or die $!;
my $counter = 0;
while(my $line = <COV>){
    chomp $line;
    (my $chr, my $pos, my $cov) = split "\t", $line;
    print $line,"\t";
    for(my $i = 0; $i <= $#{$regions{$chr}}; $i++){
	# if current position is within region add coverage
	my $end = $regions{$chr}[$i]{'end'};
	my $begin = $regions{$chr}[$i]{'begin'};
	print $begin,' ', $end,"\n";
	if($pos <= $begin){
	    last;
	}elsif($pos > $end){
	    next;
	}else{
	    $regions{$chr}[$i]{'cov'} += $cov;
	    print '   in', $regions{$chr}[$i]{'cov'},"\n";
	    last; # save time as we expect sorted
	}
	# if not continue
#use Data::Dumper; print Dumper %regions;die;
    }
    $counter++;
    last if $counter >10000;
}
close COV;
#use Data::Dumper; print Dumper %regions;die;
foreach my $chr (sort keys %regions){
    foreach my $be (@{$regions{$chr}}){
	printab($chr, $be->{'begin'}, $be->{'end'}, $be->{'coverage'}, $be->{'name'});
    }	
}
