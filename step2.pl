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

if(0){
    system "./msaLocal.pl > msaLocal.log";
    system "Rscript --vanilla ./filterProtAlign.R";
    my $cmdl = './prot2nucl.pl';
    foreach my $sp (keys %$nr){
	my %pab = %{$nr->{$sp}}; #de-reference for less typing
	$cmdl .= ' -cds '.catfile($pab{'datapath'}, $pab{'cdsF'});
    }
    $cmdl .= ' -cds /home/dsellis/data/IES/thermophila/gene/T_thermophila_June2014_CDS.fasta'; #add Tth
#$cmdl .= ' -cds /home/dsellis/data/IES/thermophila/gene/T_thermophila_June2014.gene.fa'; #add Tth
    $cmdl .= ' ~/data/IES/analysis/msas/filtered/*.aln.fa';
    run($cmdl, 0);
}

run('./runGblocks.pl ~/data/IES/analysis/msas/filtered/ > ~/data/IES/analysis/msas/filtered/gblocks.log', 1);
run('./parseGblocks.pl  ~/data/IES/analysis/msas/filtered/ > ~/data/IES/analysis/msas/filtered/parsegblocks.log', 1);

if(0){
    system "./preparePhyldog.pl";
}
