#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

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

# read iesInGenes.msa.tab
#geneFamily      gene    beginMSA        endMSA  beginGene       endGene ies
# 645     PBIA.V1_4.1.P03540030   921     922     849     850     pbi.MICA.0354.59487
# 645     PBIA.V1_4.1.P03540030   921     922     849     850     pbi.MICA.0354.59487
# 645     PBIA.V1_4.1.P03540030   1277    1278    1202    1203    pbi.MICA.0354.59064
# 645     PBIA.V1_4.1.P03540030   1277    1278    1202    1203    pbi.MICA.0354.59064

# read homColinBiB.be
#geneFamily beginMSA endMSA IES
# 10609   132     134     pca.MICA.0028.262818
# 1061    437     439     pca.MICA.0041.91286
# 1061    746     748     pso.MICA.108.107498
# 1061    982     984     pbi.MICA.0041.238017,poc.MICA.122.143635,ppr.MICA.133.236216
# 1061    1428    1430    pbi.MICA.0041.237685,poc.MICA.122.143967
# 1061    2957    2959    pso.MICA.124.259846

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
	    push @{$genesH{$gene}}, $ies;
	}else{
	    $genesH{$gene} = [$ies];
	}
}
close ING;

my $id = 0;
open HCB, $homColinBF or die $!; # only look at homologous columns that are fully within Gblocks

while(my $line = <HCB>){
	chomp $line;
	(my $geneFamily, my $beginMSA, my $endMSA, my $iesL) = split " ", $line;
	$id++; # each row is an entry
	$beginMSA++; # transform to 1-base index
#	print $id,"\t$geneFamily\t";
	my @iesList = split ",", $iesL;
#	print "@iesList\n";
	my %iesLH = map{$_ => $_ } @iesList; # list of IES in homologous column
	# for each gene in the gene family check if it has an IES in the list
##	printab($id, $geneFamily, $beginMSA, $endMSA);
	foreach my $gene (@{$gfH{$geneFamily}}){
            # homIESid geneFamily beginMSA endMSA gene beginGene(s) endGene(s) IESid
#	    print "\t", $gene,' ';
	    # does the gene have IESs
	    my $beginGene = 'NA'; # genomic coordinates of IESs of current homologous IES group
	    my $endGene   = 'NA'; # by default NA if there are no IES annotated
	    my $iesId     = 'NA';
	    if(defined($genesH{$gene})){
#		print 'has ies ';
		my @geneHomIES;
		my @iesInHomIES; # id of IES(s) within current homologous IES
		my @beginMSAInHomIES; # MSA coordinates of beginning for corresponding IES
		my @endMSAInHomIES;
		my @beginGeneInHomIES;
		my @endGeneInHomIES;
		foreach my $iesInGene (@{$genesH{$gene}}){ # for all IES in the current gene
		    if($iesLH{$iesInGene}){                # is any found in the IES list of the current homologous group?
			push @geneHomIES, $iesInGene;
			push @iesInHomIES, $iesLH{$iesInGene};
			push @beginMSAInHomIES, @{$iesH{$iesLH{$iesInGene}}{'beginMSA'}};
			push @endMSAInHomIES, @{$iesH{$iesLH{$iesInGene}}{'endMSA'}};
			push @beginGeneInHomIES, @{$iesH{$iesLH{$iesInGene}}{'beginGene'}};
			push @endGeneInHomIES, @{$iesH{$iesLH{$iesInGene}}{'endGene'}};
		    }
		}
		if(@geneHomIES){ # if gene has an IES from the current homologous IES group
		    $beginGene = join(',',@beginGeneInHomIES);
		    $endGene   = join(',', @endGeneInHomIES);
		    $iesId     = join(',', @geneHomIES);
		}
	    }
	    printab($id, $geneFamily, $beginMSA, $endMSA, $gene, $beginGene, $endGene, $iesId);
#		print 'no ies';
	    #	    print "\n";

	}
}
close HCB;
