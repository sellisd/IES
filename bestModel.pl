#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use lib'.';
use functions;
my $help;
my $nex;
my $table;
my $model;

my $usage = <<HERE;

Create a nexus partition file and a summary file of the best model for each gene family
usage bestModel.pl [OPTIONS] INPUTFILE(S)
where OPTIONS can be:
  -nex:    nexus output file
  -table:  table with gene family information
  -model:  overwrite best model and use given model for all partitions
  -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?'  => \$help,
			      'nex=s'   => \$nex,
			      'table=s' => \$table,
			      'model=s' => \$model
		   ));
die $usage if $#ARGV < 0;
die $usage if $help;

open NEX, '>', $nex or die $!;
print NEX "#nexus\n";
print NEX "begin sets;\n";

my $counter = 1;
my @partitionModels;
open TB, '>', $table or die $!;
print TB 'InputFile', "\t", 'geneFamily', "\t", 'hasTth', "\t", 'bestModel', "\n";
foreach my $inputFile (@ARGV){
    print $inputFile,"\n";
    my $hasTth = 0;
    my $in = Bio::SeqIO->new(-file   => $inputFile,
			     -format => 'Fasta');
    while(my $seqO = $in->next_seq()){
	my $species = $seqO->id;
	$hasTth = 1 if $species eq 'Tetrahymena_thermophila';
    }
    my $partitionString = 'part'.$counter;
    print NEX "\tcharset ", $partitionString, " = ", $inputFile, ":*;\n";
    my $bestModel;
    open IQ, $inputFile.'.iqtree' or die $!;
    while (my $line = <IQ>){
	if ($line =~ /^Best-fit model according to BIC:\s*(.*)\s*$/){
	    $bestModel = $1;
	}
    }
    $bestModel = $model if $model;
    $inputFile =~ /cluster\.(\d+)\.nucl\.fa\.renamed$/ or die $inputFile;
    my $geneFamily = $1;
    print TB $inputFile, "\t", $geneFamily, "\t", $hasTth, "\t", $bestModel, "\n";
    push @partitionModels, $bestModel.':'.$partitionString;
    close IQ;
    $counter++;
}
print NEX "charpartition byGene = ",join(", ", @partitionModels),";";
print NEX "end;\n";
close NEX;
close TB;
