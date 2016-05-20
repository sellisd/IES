#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# Reads homColinB.be       homologous columns bed file with information on homologous IES groups locations in MSA coordinates (ranges)
#  and  iesInGenes.msa.tab location of IES in genes in both MSA and gene coordinates, floating IES are in separate rows

my $iesInGenesF = '/home/dsellis/data/IES/analysis/tables/iesInGenes.msa.tab';
my $homColinBF = '/home/dsellis/data/IES/analysis/tables/homColinB.be';
my $gfs = '/home/dsellis/data/IES/analysis/iesdb/geneFamilydb.dat';

my %gfH;
# read genes in gene family
open GF, $gfs or die $!;
while(my $line = <GF>){
    chomp $line;
    (my $geneFamily, my $seqNo, my $avPairId, my $genesL) = split " ", $line;
    my @geneList = split ",", $genesL;
    die if defined( $gfH{$geneFamily});
    $gfH{$geneFamily} = [@geneList];
}
close GF;

open ING, $iesInGenesF or die $!;
my %iesH; # ies1=>{}, ies2=>{}     All IES within genes
my %genesH; # gene=>[ies1,ies2...] All genes with an IES
readline(ING); # header
while(my $line = <ING>){
	chomp $line;
	(my $geneFamily, my $gene, my $beginMSA, my $endMSA, my $beginGene, my $endGene, my $ies) = split " ", $line;
	if(defined($iesH{$ies})){
	    push @{$iesH{$ies}{'beginMSA'}}, $beginMSA; # if floating will have multiple locations
	    push @{$iesH{$ies}{'endMSA'}}, $endMSA;
	    push @{$iesH{$ies}{'beginGene'}}, $beginGene;
	    push @{$iesH{$ies}{'endGene'}}, $endGene;
	}else{
	    $iesH{$ies}={
		'geneFamily' => $geneFamily,
		'gene'       => $gene,
		'beginMSA'   => [$beginMSA],
		'endMSA'     => [$endMSA],
		'beginGene'  => [$beginGene],
		'endGene'    => [$endGene]
	    }
	}
	if(defined($genesH{$gene})){
	    $genesH{$gene}{$ies} = 1;
	}else{
	    $genesH{$gene} = {$ies => 1}; # if floating a hash will keep a singe entry for each IES
	}
}
close ING;

my $id = 0;
open HCB, $homColinBF or die $!; # only look at homologous columns that are fully within Gblocks
printab('id', 'geneFamily', 'beginMSArange', 'endMSArange', 'gene', 'beginGene', 'endGene', 'beginMSA', 'endMSA', 'ies'); # gene has other IESs, not in current homologous IES group
while(my $line = <HCB>){
	chomp $line;
	(my $geneFamily, my $beginMSA, my $endMSA, my $iesL) = split " ", $line;
	$id++; # each row is an entry
	$beginMSA++; # transform to 1-base index
	my @iesList = split ",", $iesL;
	my %iesLH = map{$_ => $_ } @iesList; # list of IES in homologous column
	# for each gene in the gene family check if it has an IES in the list
	foreach my $gene (@{$gfH{$geneFamily}}){
	    # does the gene have IESs
	    my $beginGene = 'NA'; # genomic coordinates of IESs of current homologous IES group
	    my $endGene   = 'NA'; # by default NA if there are no IES annotated
	    my $iesId     = 'NA';
	    if(defined($genesH{$gene})){ # has ies
		my $found = 0;
		foreach my $iesInGene (keys %{$genesH{$gene}}){ # for all IES in the current gene
		    if($iesLH{$iesInGene}){                # is any found in the IES list of the current homologous group?
			$iesId = $iesLH{$iesInGene};
			my $locations = $#{$iesH{$iesId}{'beginGene'}} + 1;
			for(my $c = 0; $c < $locations; $c++){
			    printab($id, $geneFamily, $beginMSA, $endMSA, $gene, ${$iesH{$iesId}{'beginGene'}}[$c], ${$iesH{$iesId}{'endGene'}}[$c], ${$iesH{$iesId}{'beginMSA'}}[$c], ${$iesH{$iesId}{'endMSA'}}[$c],$iesId);
			}
			$found = 1;
		    }
		}
		if($found == 0){
		    printab($id, $geneFamily, $beginMSA, $endMSA, $gene, 'NA', 'NA', 'NA', 'NA', 'NA'); # gene has other IESs, not in current homologous IES group
		}
	    }else{
		printab($id, $geneFamily, $beginMSA, $endMSA, $gene, 'NA', 'NA', 'NA', 'NA', 'NA'); # gene does not have IESs
	    }
	}
}
close HCB;
