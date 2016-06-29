#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;

my $gfF = '/home/dsellis/data/IES/analysis/iesdb/geneFamilydb.dat';
my $iqtreeB = '/home/dsellis/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 '; # binary location
my $pathIN = '/home/dsellis/data/IES/analysis/msas/filtered/';
my $pathOUT = '/home/dsellis/data/IES/analysis/sgf/';
my $nex =  '/home/dsellis/data/IES/analysis/sgf/part.nex';
my $concatF = catfile($pathOUT, 'concat.fa');
make_path($pathOUT) unless -e $pathOUT;

my @selectedGroups;
my @gtF; # gene tree files selected
my @partitionModels;

open IN, $gfF or die $!;
readline(IN); #header
# find which gene families have only one memeber of each species
while (my $line = <IN>){
    chomp $line;
    (my $id, my $seqNo, my $avPairId, my $genes, my $pprGenes, my $pbiGenes, my $pteGenes, my $ppeGenes, my $pseGenes, my $pocGenes, my $ptrGenes, my $psoGenes, my $pcaGenes) = split " ", $line;
    if($pprGenes == 1 and $pbiGenes == 1 and  $pteGenes == 1 and $ppeGenes == 1 and $pseGenes == 1 and $pocGenes == 1 and $ptrGenes == 1 and $psoGenes == 1 and $pcaGenes == 1){
	my $fileName = 'cluster.'.$id.'.nucl.fa';
	push @selectedGroups, $id;
	push @gtF, catfile($pathIN, $fileName);
    }
}
close IN;

# concatenate alignment files
my $cmdl = './concat.pl '.join(' ', @gtF);
run($cmdl, 1);


$cmdl = $iqtreeB.' -s '.$concatF.
    ' -m GTR+G{1.0}'.
    ' -redo'.
    ' -b 1000';
run($cmdl, 0);
