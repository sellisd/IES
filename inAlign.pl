#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Bio::AlignIO;
use File::Spec::Functions qw(catfile);

# read gene - gene family

my $silixOutF = '/home/dsellis/data/IES/analysis/allvsall/blastout/silix.output';
my $iesig = '/home/dsellis/data/IES/analysis/tables/ppe.iesInGenes';
my $alnPath = '/home/dsellis/data/IES/analysis/msas/filtered/';

my $prefixesR = initF();

open S, $silixOutF or die $!;
my %gf; #gene families
while(my $line = <S>){
    (my $geneFamily, my $protName) = split " ", $line;
    my $geneName = prot2gene($protName, $prefixesR);
    $gf{$geneName} = $geneFamily;
}
close S;

# read ies in gene
open IG, $iesig or die $!;
while(my $line = <IG>){
    (my $gene, my $start, my $end, my $ies) = split " ", $line;
    my $geneFamily = $gf{$gene};
# save data structure per gene family
    die;
}
close IG;

#once all are read go trhough alignments one by one
    my $alnF = catfile($alnPath, 'cluster.'.$geneFamily.'.nucl.fa');
    next unless (-e $alnF); # skip filtered files
    my $alnIO = Bio::AlignIO->new(-file => $alnF,
				  -format => 'fasta');
    while(my $alnO = $alnIO->next_aln){ # open alignment and find position of ies in alignment
	my $seqNo = $alnO->num_sequences;
	print $alnF, ' ', $seqNo, ' ';
	print $alnO->percentage_identity(), "\n";
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $id = $seq->id();

    }
