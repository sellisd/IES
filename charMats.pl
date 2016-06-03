#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use lib'.';
use functions;

# load gene/protein notation
my $prefixesR = initF();

# prepare character matrices

# read list of genes which have too low read coverage in the Mic genome and thus any IES detected in them or homologous IES present in them is of low confidence
my $path = '/home/dsellis/data/IES/analysis/tables/';
my @gF = qw/gpbi.filt gpca.filt gpoc.filt gppe.filt gppr.filt gpse.filt gpso.filt gpte.filt gptr.filt/; # files with filtered genes 

#my $covgeneF = '/home/dsellis/data/IES/analysis/tables/pca.gene.cov';
my %h; # list of filtered genes
foreach my $file (@gF){
    open CV, catfile($path, $file) or die $!;
    while(my $gene = <CV>){
	chomp $gene;
	die if defined ($h{$gene}); #they should be unique records
	$h{$gene} = 1;
    }
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
    # if gene in list print ? and keep track of whether it was present or absent
    # else print 1 if present 0 if absent
    if(defined($h{$gene})){
	$ies = '?';
    }else{
	if($ies eq 'NA'){
	    $ies = 0;
	}else{
	    $ies = 1;
	}
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

