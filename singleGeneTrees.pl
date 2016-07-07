#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use File::Copy;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new(7);

my $notationF = '/home/dsellis/data/IES/analysis/notation.csv';
my $gfF = '/home/dsellis/data/IES/analysis/iesdb/geneFamilydb.dat';
my $iqtreeBP = '/home/dsellis/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 '; # binary location
my $iqtreeB =  '/home/dsellis/tools/iqtree-1.4.2-Linux/bin/iqtree';
my $pathIN = '/home/dsellis/data/IES/analysis/msas/filtered/';
my $pathOUT = '/home/dsellis/data/IES/analysis/sgf/';
my $nex =  '/home/dsellis/data/IES/analysis/sgf/part.nex';
my $concatF = catfile($pathOUT, 'concat.fa');
my $bmF = '/home/dsellis/data/IES/analysis/tables/bestModels.tab'; # best models for each single gene family
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
$cdsF .= ' -cds /home/dsellis/data/IES/genomicData/thermophila/gene/T_thermophila_June2014_CDS.fasta'; #add Tth

run('./prot2nucl.pl -noterm'.$cdsF.' ~/data/IES/analysis/sgf/*.aln.fa', 1);

# rename sequences
run('./nameReplaceAlign.pl ~/data/IES/analysis/sgf/cluster.*.nucl.fa', 1);

# find best model for each gene family and infer gene tree
opendir DH, $pathOUT or die $!;
my @nuclAlnF = grep {/.*\.nucl\.fa\.renamed$/} readdir(DH);
close DH;
foreach my $file (@nuclAlnF){
    my $pid = $pm->start and next;
    my $cmdl = "$iqtreeB -s ".catfile($pathOUT, $file).
	' -st CODON6 -bb 1000'.
	' -m TESTNEW';
    run($cmdl, 1);
    $pm->finish;
}
$pm->wait_all_children;

# build table with best models for each partition
run('./bestModel.pl -nex ~/data/IES/analysis/sgf/concat.nexus -table '.$bmF.' ~/data/IES/analysis/sgf/cluster.*.nucl.fa.renamed', 1);

# infer concatenated (species) tree with partitions and -testmerge
run($iqtreeBP.' -bb 1000 -st CODON6 -m TESTNEWMERGE -spp  ~/data/IES/analysis/sgf/concat.nexus > ~/data/IES/analysis/log/concat.log', 0);

# infer concatenated (species) tree with simple model
run('./bestModel.pl -model GTR+G{1.0} -nex ~/data/IES/analysis/sgf/concatSimple.nexus -table '.$bmF.' ~/data/IES/analysis/sgf/cluster.*.nucl.fa.renamed', 0);
run($iqtreeBP.' -bb 1000 -spp  ~/data/IES/analysis/sgf/concatSimple.nexus > ~/data/IES/analysis/log/concatSimple.log', 0);

# compare gene trees with simple model and partitioned species tree topology

