#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Getopt::Long;
use Bio::Tools::GFF;

my $help;

my $silixOutF = '/home/dsellis/data/IES/analysis/allvsall/blastout/silix.output';
my $gff = '/home/dsellis/data/IES/primaurelia/gene/pprimaurelia_AZ9-3_annotation_v1.0.gff3';
my $usage = <<HERE;

combine the result of analyses to create the IES table for iesDB
usage genedbTable.pl [OPTIONS]
where OPTIONS can be:
    -silixout : file with the silix output containing the gene family assignments
    -gff :      gff file with gene annotation
    -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?'   => \$help,
			      'silixout=s'  => \$silixOutF,
			      'gff=s'  => \$gff,
		   ));
die $usage if $help;

# read silix output and get genes per cluster
open S, $silixOutF or die $!;
my %gf; #gene families
my $prefixesR = initF();

while(my $line = <S>){
    (my $geneFamily, my $protName) = split " ", $line;
    my $geneName = prot2gene($protName, $prefixesR);
    $gf{$geneName} = $geneFamily;
}
close S;

# read gff file and build CDS in each gene
my $gffI = Bio::Tools::GFF->new('-file' => $gff,
 				'-gff_version' => 3);
my %genes;
printab('id', 'scaffold', 'geneName', 'strand', 'start', 'end', 'geneFamily');
while(my $feature = $gffI->next_feature()){
    my $scaffold = $feature->seq_id;
    my $start = $feature->start;
    my $end = $feature->end;
    my $id = $feature->primary_id;
    my $strand = $feature->strand;
    if($feature->primary_tag() eq 'CDS'){
	my @parent = $feature->get_tag_values('Parent');
	my $geneName = X2Y($parent[0], $prefixesR, 'T', 'G'); #gene from transcript
	my $geneFamily = $gf{$geneName};
	if(!defined($geneFamily)){
	    # e.g in small scaffold
	    $geneFamily = 'NA';
	}
	printab($id, $scaffold, $geneName, $strand, $start, $end, $geneFamily);
    }
}
