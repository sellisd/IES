#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# combine information to create the gene family table
my $silixOutF = '/home/dsellis/data/IES/analysis/allvsall/blastout/silix.output';
my $protAlignStatsF = '/home/dsellis/data/IES/analysis/msas/protAlignStats.dat';

my $outF = '/home/dsellis/data/IES/analysis/iesdb/geneFamilydb.dat';

#read silix output
print "reading SiLiX gene families\n";

my %geneFamily;
my $prefixesR = initF();
open S, $silixOutF or die $!;
while(my $line = <S>){
    (my $gf, my $protName) = split " ", $line;
    my $geneName = prot2gene($protName, $prefixesR);
    if(defined($geneFamily{$gf})){
	push @{$geneFamily{$gf}{'genes'}}, $geneName;
    }else{
	$geneFamily{$gf}{'genes'} = [$geneName];
    }
}
close S;

# read filtering stats
print "reading protein alignment statistics\n";

open FS, $protAlignStatsF or die $!;
readline(FS); #header
while(my $line = <FS>){
    chomp $line;
    (my $id, my $seqNo, my $avpid) = split " ", $line;
    $id =~ s/^cluster\.(\d+)\.aln\.fa$/$1/ or die $line;
    $geneFamily{$id}{'seqNo'}  = $seqNo;
    $geneFamily{$id}{'avPairId'} = $avpid;
}
close FS;

open OUT, '>', $outF or die $!;
print OUT join("\t", ("id", "seqNo", "avPairId", "genes", "pprGenes", "pbiGenes", "pteGenes", "ppeGenes", "pseGenes", "pocGenes", "ptrGenes", "psoGenes", "pcaGenes")), "\n";
foreach my $gf (sort keys %geneFamily){
    print OUT $gf,"\t";
    if(defined($geneFamily{$gf}{'seqNo'})){ # not filtered !
	print OUT $geneFamily{$gf}{'seqNo'},"\t";
	print OUT $geneFamily{$gf}{'avPairId'},"\t";
    }else{
	print OUT 'NA',"\t";
	print OUT 'NA',"\t";
    }
    print OUT join(',', @{$geneFamily{$gf}{'genes'}}), "\t";
    my $histR = tabGenes($geneFamily{$gf}{'genes'});
    print OUT $histR->{'Paramecium_primaurelia'}, "\t";
    print OUT $histR->{'Paramecium_biaurelia'}, "\t";
    print OUT $histR->{'Paramecium_tetraurelia'}, "\t";
    print OUT $histR->{'Paramecium_pentaurelia'}, "\t";
    print OUT $histR->{'Paramecium_sexaurelia'}, "\t";
    print OUT $histR->{'Paramecium_octaurelia'}, "\t";
    print OUT $histR->{'Paramecium_tredecaurelia'}, "\t";
    print OUT $histR->{'Paramecium_sonneborni'}, "\t";
    print OUT $histR->{'Paramecium_caudatum'}, "\t";
    print OUT "\n";
}
close OUT;

# #output
# id seqNo avPairId Genes geneNumber
sub tabGenes{
# tabulate genes by species
    my $genesR = shift @_;
    my %hist = (
	'Paramecium_primaurelia'   => 0,
	'Paramecium_biaurelia'     => 0,
	'Paramecium_tetraurelia'   => 0,
	'Paramecium_pentaurelia'   => 0,
	'Paramecium_sexaurelia'    => 0,
	'Paramecium_octaurelia'    => 0,
	'Paramecium_tredecaurelia' => 0,
	'Paramecium_sonneborni'    => 0,
	'Paramecium_caudatum'      => 0
	);
    foreach my $gene (@$genesR){
	$hist{gene2species($gene)}++;
    }
    return \%hist;
}
