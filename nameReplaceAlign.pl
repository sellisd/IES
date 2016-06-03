#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use lib'.';
use functions;
use Getopt::Long;
my $help;
my $nex;
my $usage = <<HERE;

Replace gene names with species names in multiple alignment files and create a nexus partition file at the same time
usage nameReplaceAlign.pl [OPTIONS] INPUTFILE(S)
where OPTIONS can be:
  -nex:    nexus output file
  -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'nex=s'  => \$nex
		   ));
die $usage if $#ARGV < 0;
die $usage if $help;

my @species = qw/
Paramecium_primaurelia
Paramecium_biaurelia
Paramecium_tetraurelia
Paramecium_pentaurelia
Paramecium_sexaurelia
Paramecium_octaurelia
Paramecium_tredecaurelia
Paramecium_sonneborni
Paramecium_caudatum
Tetrahymena_thermophila
/;

open NEX, '>', $nex or die $!;
print NEX "#nexus\n";
print NEX "begin sets;\n";
my $counter = 1;
my @partitionModels;
foreach my $inputFile (@ARGV){
    print $inputFile,"\n";
    my $in = Bio::SeqIO->new(-file   => $inputFile,
			     -format => 'Fasta');
    my $outputFile = $inputFile.'.renamed';
    my $partitionString = 'part'.$counter;
    print NEX "\tcharset ", $partitionString, " = ", $outputFile, ":*;\n";
    my $out = Bio::SeqIO->new(-file => '>'.$outputFile,
			      -format => 'Fasta');
    while(my $seqO = $in->next_seq()){
	my $species = gene2species($seqO->id);
	my $newSeqO = Bio::Seq->new(-id => $species,
				    -seq => $seqO->seq());
	$out->write_seq($newSeqO);
    }
    $counter++;
    # find best model
    my $bestModel;
    open IQ, $inputFile.'.iqtree' or die $!;
    while (my $line = <IQ>){
	if ($line =~ /^Best-fit model according to BIC:\s*(.*)\s*$/){
	    $bestModel = $1;
	}
    }
    push @partitionModels, $bestModel.':'.$partitionString;
    close IQ;
}
print NEX "charpartition byGene = ",join(", ", @partitionModels),";";
print NEX "end;\n";
close NEX;
