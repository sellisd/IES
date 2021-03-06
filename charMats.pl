#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile catdir);
use lib'.';
use functions;

my $opt = loadUserOptions;
my $basePath   = $$opt{'basePath'};
my $tablesP    = catdir($basePath, 'analysis', 'tables');

# load gene/protein notation
my $prefixesR = initF();

# prepare character matrices

# read list of genes which have sufficient read coverage in the Mic genome and thus any IES detected in them or homologous IES present in them is of confidence
my @gF = qw/gpbi.filt gpca.filt gpoc.filt gppe.filt gppr.filt gpse.filt gpso.filt gpte.filt gptr.filt/; # files with filtered genes 

#my $covgeneF = '/home/dsellis/data/IES/analysis/tables/pca.gene.cov';
my %h; # list of filtered genes
foreach my $file (@gF){
    open CV, catfile($tablesP, $file) or die $!;
    while(my $gene = <CV>){
	chomp $gene;
	die if defined ($h{$gene}); #they should be unique records
	$h{$gene} = 1;
    }
}
close CV;

#floatingIES ptr.MICA.7180000129095.207151

my $homIESdb = catfile($tablesP, 'homIES.tab');
open IN, $homIESdb or die "$! $homIESdb";
readline(IN); # header
my %charMats;
while(my $line = <IN>){
    chomp $line;
    (my $id, my $geneFamily, my $beginMSArange, my $endMSArange, my $gene, my $beginGene, my $endGene, my $beginMSA, my $endMSA, my $iesId) = split " ", $line;
    # if T. thermophila print 0 (absent)
    # if gene is not in list print ?
    # else print 1 if present 0 if absent
    my $ies; # presence absence or unknown
    if(gene2species($gene) eq 'Tetrahymena_thermophila'){
	$ies = 0;
    }else{
	if(defined($h{$gene})){
	    if($iesId eq 'NA'){
		$ies = 0;
	    }else{
		$ies = 1;
	    }
	}else{
	    $ies = '?';
	}
    }
    $charMats{$geneFamily}{$id}{$gene} = {'begin'    => $beginMSArange,
					  'end'      => $endMSArange,
					  'ies'      => $ies,
					  'iesId'    => $iesId,
					  'beginMSA' => $beginMSA,
					  'endMSA'   => $endMSA};
}
close IN;

# character matrix output, translate geneId to protId to keep it comparable with the trees
printab(('cluster', 'column', 'geneId', 'begin', 'end', 'ies', 'iesId', 'beginMSA', 'endMSA'));
foreach my $geneFamily (sort {$a<=>$b} keys %charMats){
    foreach my $id (sort {$a <=> $b} keys %{$charMats{$geneFamily}}){
	foreach my $gene (sort keys  %{$charMats{$geneFamily}{$id}}){
	    printab($geneFamily, $id, X2Y($gene, $prefixesR, 'G', 'P'), $charMats{$geneFamily}{$id}{$gene}{'begin'},
		    $charMats{$geneFamily}{$id}{$gene}{'end'},
		    $charMats{$geneFamily}{$id}{$gene}{'ies'},
		    $charMats{$geneFamily}{$id}{$gene}{'iesId'},
		    $charMats{$geneFamily}{$id}{$gene}{'beginMSA'},
		    $charMats{$geneFamily}{$id}{$gene}{'endMSA'}
		);
	}
    }
}

