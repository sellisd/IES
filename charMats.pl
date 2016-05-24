#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# load gene/protein notation
my $prefixesR = initF();

# prepare character matrices

# read read coverage per gene (Only for P. caudatum) if a gene has on average coverage more than cutoff (15x) we are confident on the absence of IES, otherwise it should be coded as ?

my $covgeneF = '/home/dsellis/data/IES/analysis/tables/pca.gene.cov';
my %h;
open CV, $covgeneF or die $!;
while(my $line = <CV>){
    chomp $line;
    (my $gene, my $coverage) = split " ", $line;
    $h{$gene} = $coverage;
}
close CV;

#floatingIES ptr.MICA.7180000129095.207151

my $homIESdb = '/home/dsellis/data/IES/analysis/iesdb/homIESdb.dat';
open IN, $homIESdb or die $!;
readline(IN); # header
my %charMats;
while(my $line = <IN>){
    chomp $line;
    (my $id, my $geneFamily, my $beginMSArange, my $endMSArange, my $gene, my $beginGene, my $endGene, my $beginMSA, my $endMSA, my $ies) = split " ", $line;
    if($ies eq 'NA'){
	if(defined($h{$gene})){ # if no coverage data then this is a species with no low-coverage problem
	    if($h{$gene} < 15){
		$ies = '?';
	    }elsif($h{$gene} >= 15){
		$ies = 0;
	    }else{
		die;
	    }
	}else{
	    $ies = 0;
	}
    }else{
	$ies = 1;
    }
    $charMats{$geneFamily}{$id}{$gene} = {'begin' => $beginMSArange,
					  'end'   => $endMSArange,
					  'ies'   => $ies};
}
close IN;

# character matrix output, translate geneId to protId to keep it comparable with the trees
printab(('cluster', 'column', 'geneId', 'begin', 'end', 'ies'));
foreach my $geneFamily (sort {$a<=>$b} keys %charMats){
    foreach my $id (sort {$a <=> $b} keys %{$charMats{$geneFamily}}){
	foreach my $gene (sort keys  %{$charMats{$geneFamily}{$id}}){
	    printab($geneFamily, $id, X2Y($gene, $prefixesR, 'G', 'P'), $charMats{$geneFamily}{$id}{$gene}{'begin'},
		    $charMats{$geneFamily}{$id}{$gene}{'end'},
		    $charMats{$geneFamily}{$id}{$gene}{'ies'});
	}
    }
}

