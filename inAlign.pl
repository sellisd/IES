#!/usr/bin/perl
use warnings;
use strict;

# read gene - gene family
# read ies in gene
# open alignment and find position of ies in alignment

my %geneH;
my %cdsH;
my $genedb = '/home/dsellis/data/IES/analysis/iesdb/ppr.cdsdb';
my $iesinF = '/home/dsellis/data/IES/analysis/bed/ppr.IESin.be';

open GB, $genedb or die $!;
readline(GB);
while (my $line = <GB>){
    chomp $line;
    (my $id, my $scaffold, my $geneName, my $strand, my $start, my $end, my $geneFamily) = split " ", $line;
    $cdsH{$id} = $geneName;
    if(defined($geneH{$geneName})){
	push @{$geneH{$geneName}{'cdsStart'}}, $start;
	push @{$geneH{$geneName}{'cdsEnd'}}, $end;
	push @{$geneH{$geneName}{'cds'}}, $id;
    }else{
	$geneH{$geneName} = {'scaffold'   => $scaffold,
			     'strand'     => $strand,
			     'geneFamily' => $geneFamily,
			     'cdsStart'   => [$start],
			     'cdsEnd'     => [$end], 
			     'cds'        => [$id]};
    }
}
close GB;

open IN, $iesinF or die $!;
while (my $line = <IN>){
    chomp $line;
    (my $scaffold, my $start, my $end, my $name, my $isFloating, my $inType, my $inScaffold, my $inStart, my $inEnd, my $name1, my $overlap) = split " ", $line;
    next unless $inType eq 'cds';
    my $gene = $cdsH{$name1}; #gene in which ies is found
    my $strand = $geneH{$gene}{'strand'};
    my @cdsNames = $geneH{$gene}{'cds'};
    my @cdsStart = $geneH{$gene}{'cdsStart'};
    my @cdsEnd = $geneH{$gene}{'cdsEnd'};
    die unless $scaffold eq $geneH{$geneName}{'scaffold'}; #sanity check
    my $cumLength = 0;
    my $startTC;
    my $endTC;
    my $prevCdsStart = 0; # sanity check, cds coordinates should be sorted
    if($strand == '1'){
#	print "ies $name is in $name1 which is in $gene\n";
	for(my $i = 0; $i <= $#cdsNames; $i++){
	    die unless ($cdsStart[$i] > $prevCdsStart);
	    $prevCdsStart = $cdsStart[$i];
	    if($name1 eq $cdsNames[$i]){
# add partial length
		$startTC = $cumLength + $start - $inStart; # ies start in transcript coordinates zero-based
		$endTC = $cumLength + $end - $inStart;
		last;
	    }else{
		$cumLength += $cdsEnd[$i] - $cdsStart[$i] + 1; # coordinates from the gff file
	    }
	}
    }elsif($strand == '-1'){
	for(my $i = $#cdsNames; $i >= 0; $i--){ #loop backwards
	    die unless ($cdsStart[$i] < $prevCdsStart);
	    $prevCdsStart = $cdsStart[$i];
	    if($name eq $cdsNames[$i]){
		$startTC = $cumLength + $inEnd - $start;
		$endTC = $cumLength + $inEnd - $end;
	    }else{
		$cumLength += $cdsEnd[$i] - $cdsStart[$i] + 1;
	    }	    
	}
    }else{
	die;
    }
if in strand is forward
add CDS length up to current cds and then add partial
if -
add remaining
}
   

close IN;
# # read IESin.be
# ies in cds
# find in which gene
# find gene orientation & gene family
# calculate distance from begin/end of cds based on overlap
# add total length of remaining cds of gene:
# sort cds by begin find all following/preceding add length

# for each ies location find in which CDS it is (if in any)
#   then find in which gene this CDS is in
  


# read alignments
# find location of ies in alignments
