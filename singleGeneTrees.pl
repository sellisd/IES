#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use File::Copy;

my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};

my $notationF = catfile($basePath, 'analysis', 'notation.csv');
my $gfF       = catfile($basePath, 'analysis', 'iesdb', '/geneFamilydb.dat');
my $iqtreeBP  = '/home/dsellis/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 '; # binary location
my $pathIN    = catfile($basePath, 'analysis', 'msas', 'filtered');
my $pathOUT   = catfile($basePath, 'analysis', 'sgf');
my $nex       = catfile($basePath, 'analysis', 'sgf', 'part.nex');
my $bmF       = catfile($basePath, 'analysis', 'tables', 'bestModels.tab'); # best models for each single gene family
my $tthF      = catfile($basePath, 'genomicData', 'thermophila' ,'gene', 'T_thermophila_June2014_CDS.fasta')
my $sgfP      = catfile($basePath, 'analysis', 'sgf')
my $concatF   = catfile($pathOUT, 'concat.fa');

make_path($pathOUT) unless -e $pathOUT;

my @selectedGroups;
my @gtF; # gene tree files selected
my @partitionModels;

open IN, $gfF or die $!;
readline(IN); #header
# find which gene families have only one memeber of each species

while (my $line = <IN>){
    chomp $line;
    (my $id, my $seqNo, my $avPairId, my $genes, my $pprGenes, my $pbiGenes, my $pteGenes, my $ppeGenes, my $pseGenes, my $pocGenes, my $ptrGenes, my $psoGenes, my $pcaGenes, my $tthGenes) = split " ", $line;
    if($pprGenes == 1 and $pbiGenes == 1 and  $pteGenes == 1 and $ppeGenes == 1 and $pseGenes == 1 and $pocGenes == 1 and $ptrGenes == 1 and $psoGenes == 1 and $pcaGenes == 1 and $tthGenes <=1){
	     my $fileName = 'cluster.'.$id.'.aln.fa';
	      push @selectedGroups, $id;
	       push @gtF, catfile($pathIN, $fileName);
    }
}
close IN;

# copy aminoacid alignments to sgf folder
foreach my $sgf (@selectedGroups){
    my $from = catfile($pathIN, 'cluster.'.$sgf.'.aln.fa');
    copy($from, $pathOUT);
}

# "back-translate" to nucleotide sequence without termination codon
my $nr = getNotation($notationF);
my $cdsF;
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    $cdsF .= ' -cds '.catfile($pab{'datapath'}, $pab{'cdsF'});
}
$cdsF .= ' -cds '.$tthF; #add Tth

run('./prot2nucl.pl -noterm'.$cdsF.' '.$sgfP.'/*.aln.fa', 1);

# rename sequences
run('./nameReplaceAlign.pl '.$sgfP.'/cluster.*.nucl.fa', 1);

# infer concatenated (species) tree with simple model, do not ovewrite $bmF
my $concatSimpleF = catfile($sgfP, 'concatSimple.nexus')
my $concatSimpleL = catfile($basePath, 'analysis', 'log', 'concatSimple.log')
run('./bestModel.pl -model GTR+G{1.0} -nex '.$concatSimpleF.' -table /dev/null '.$sgfP.'/cluster.*.nucl.fa.renamed', 1);
run($iqtreeBP.' -bb 1000 -spp  '.$concatSimpleF.' > '.$concatSimpleL, 1);

#move to cluster
