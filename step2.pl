#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;

my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};
my $filtscafP  = catdir($basePath, 'analysis', 'filtscaf');
my $iesdbP     = catdir($basePath, 'analysis', 'iesdb');
my $logP       = catdir($basePath, 'analysis', 'log');
my $notationF =  catfile($basePath, 'analysis', 'notation.csv');
my $tablesP    = catdir($basePath, 'analysis', 'tables');

my $nr = getNotation($notationF);


foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $cmdl = './genedb.pl'.
	' -silixout '.catfile($basePath, 'analysis', 'allvsall', 'blastout', 'silix.output').
	' -origff '.catfile($pab{'datapath'}, $pab{'geneGff'}).
	' -gff '.catfile($filtscafP, $pab{'abr'}.'.gff').
	' > '.catfile($iesdbP, $pab{'abr'}.'.genedb');
    run($cmdl, 1);
}

run('./msaLocal.pl > '.catfile($logP, 'msaLocal.log'), 1);
run("Rscript --vanilla ./filterProtAlign.R", 1);
my $cmdl;
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    $cmdl .= ' -cds '.catfile($pab{'datapath'}, $pab{'cdsF'});
}
$cmdl .= ' -cds '.catfile($basePath,'genomicData', 'thermophila', 'gene', 'T_thermophila_June2014_CDS.fasta'); #add Tth

run('./prot2nucl.pl '.$cmdl.' '.catfile($basePath, 'analysis', 'msas', 'filtered/*.aln.fa'), 1);
#run('./prot2nucl.pl -noterm'.$cmdl.' ~/data/IES/analysis/singleGene/*.aln.fa', 1);

run('./paraGblocks.pl '.catdir($basePath, 'analysis', 'msas', 'filtered').' > '.catfile($logP, 'gblocks.log'), 1);

$cmdl = './inAlign.pl'.
		' -iesig '.$tablesP.
		' -alnPath '.catfile($basePath, 'analysis', 'msas', 'filtered').
		' -out '.catfile($tablesP, 'iesInGenes.msa');
run($cmdl, 1);

$cmdl = 'sort -k 1,1 -k 2,2n '.catfile($tablesP, 'iesInGenes.msa.be').'  > '.catfile($tablesP, 'iesInGenes.msa.sorted.be');
run($cmdl, 1);
$cmdl = 'bedtools merge -i '.catfile($tablesP, 'iesInGenes.msa.sorted.be').' -c 4 -o collapse > '.catfile($tablesP, 'homColumns.be');
run($cmdl, 1);

$cmdl= './withinGblocks.pl';
run($cmdl, 1);

run('./homIESdb.pl > '.catfile($tablesP, 'homIES.tab'), 1);

$cmdl = './charMats.pl > '.catfile($iesdbP, 'charMats.tab');
run($cmdl, 1);

$cmdl = './geneFamilydb.pl';
run($cmdl, 1);

run("./preparePhyldog.pl", 1);
