#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;
my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.tab');

my $nr = getNotation($notationF);

run("./msaLocal.pl > /home/dsellis/data/IES/analysis/log/msaLocal.log", 1);
run("Rscript --vanilla ./filterProtAlign.R", 1);
my $cmdl = './prot2nucl.pl';
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    $cmdl .= ' -cds '.catfile($pab{'datapath'}, $pab{'cdsF'});
}
$cmdl .= ' -cds /home/dsellis/data/IES/genomicData/thermophila/gene/T_thermophila_June2014_CDS.fasta'; #add Tth
$cmdl .= ' ~/data/IES/analysis/msas/filtered/*.aln.fa';
run($cmdl, 1);


run('./paraGblocks.pl ~/data/IES/analysis/msas/filtered/ > ~/data/IES/analysis/log/gblocks.log', 1);

$cmdl = './inAlign.pl -iesig /home/dsellis/data/IES/analysis/tables/ -alnPath /home/dsellis/data/IES/analysis/msas/filtered/ -out /home/dsellis/data/IES/analysis/tables/iesInGenes.msa';
run($cmdl, 1);

$cmdl = 'sort -k 1,1 -k 2,2n /home/dsellis/data/IES/analysis/tables/iesInGenes.msa.be  > /home/dsellis/data/IES/analysis/tables/iesInGenes.msa.sorted.be';
run($cmdl, 1);
$cmdl = 'bedtools merge -i /home/dsellis/data/IES/analysis/tables/iesInGenes.msa.sorted.be -c 4 -o collapse > /home/dsellis/data/IES/analysis/tables/homColumns.be';
run($cmdl, 1);

$cmdl= './withinGblocks.pl';
run($cmdl, 1);

$cmdl = './geneFamilydb.pl';
run($cmdl, 1);

#'iqtree -s '..' -st CODON6 -m TESTNEWONLY'
run("./preparePhyldog.pl", 1);
