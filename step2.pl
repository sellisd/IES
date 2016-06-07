#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new(7);

my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.csv');

my $nr = getNotation($notationF);


foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $cmdl = './genedb.pl'.
	' -silixout '.'/home/dsellis/data/IES/analysis/allvsall/blastout/silix.output'.
	' -origff '.catfile($pab{'datapath'}, $pab{'geneGff'}).
	' -gff '.catfile('/home/dsellis/data/IES/analysis/filtscaf/', $pab{'abr'}.'.gff').
	' > /home/dsellis/data/IES/analysis/iesdb/'.$pab{'abr'}.'.genedb';
    run($cmdl, 1);
}

run("./msaLocal.pl > /home/dsellis/data/IES/analysis/log/msaLocal.log", 1);
run("Rscript --vanilla ./filterProtAlign.R", 1);
my $cmdl = './prot2nucl.pl';
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    $cmdl .= ' -cds '.catfile($pab{'datapath'}, $pab{'cdsF'});
}
$cmdl .= ' -cds /home/dsellis/data/IES/genomicData/thermophila/gene/T_thermophila_June2014_CDS.fasta'; #add Tth

run($cmdl.' ~/data/IES/analysis/msas/filtered/*.aln.fa', 1);
run($cmdl.' ~/data/IES/analysis/singleGene/*.aln.fa', 1);

run('./paraGblocks.pl ~/data/IES/analysis/msas/filtered/ > ~/data/IES/analysis/log/gblocks.log', 1);

$cmdl = './inAlign.pl -iesig /home/dsellis/data/IES/analysis/tables/ -alnPath /home/dsellis/data/IES/analysis/msas/filtered/ -out /home/dsellis/data/IES/analysis/tables/iesInGenes.msa';
run($cmdl, 1);

$cmdl = 'sort -k 1,1 -k 2,2n /home/dsellis/data/IES/analysis/tables/iesInGenes.msa.be  > /home/dsellis/data/IES/analysis/tables/iesInGenes.msa.sorted.be';
run($cmdl, 1);
$cmdl = 'bedtools merge -i /home/dsellis/data/IES/analysis/tables/iesInGenes.msa.sorted.be -c 4 -o collapse > /home/dsellis/data/IES/analysis/tables/homColumns.be';
run($cmdl, 1);

$cmdl= './withinGblocks.pl';
run($cmdl, 1);

$cmdl = './charMats.pl > /home/dsellis/data/IES/analysis/iesdb/charMats.tab';
run($cmdl, 1);

$cmdl = './geneFamilydb.pl';
run($cmdl, 1);

run('./homIESdb.pl > /home/dsellis/data/IES/analysis/table/homIES.dat', 1);

run("./preparePhyldog.pl", 1);

# infer single gene families phylogeny
my $singleGeneP = '/home/dsellis/data/IES/analysis/singleGene';
# rename genes in gene families
run('./nameReplaceAlign.pl ~/data/IES/analysis/singleGene/cluster.*.nucl.fa', 1);

opendir DH, $singleGeneP or die $!;
my @nuclAlnF = grep {/.*\.nucl\.fa.renamed$/} readdir(DH);
close DH;
my $iqtreeB = '/home/dsellis/tools/iqtree-1.4.2-Linux/bin/iqtree'; #binary
if(0){
foreach my $file (@nuclAlnF){
    my $pid = $pm->start and next;
#strip stop codons from alignments
# read alignment file
# for each sequence 
    my $cmdl = "$iqtreeB -s ".catfile($singleGeneP, $file).
#	' -m TESTNEWONLY -b 100';
	' -m TESTNEWONLY -redo';
    run($cmdl, 1);
    $pm->finish;
}
$pm->wait_all_children;
}
run('./bestModel.pl -nex ~/data/IES/analysis/singleGene/part.nexus -table ~/data/IES/analysis/tables/bestModels.tab ~/data/IES/analysis/singleGene/cluster.*.nucl.fa.renamed', 1);

run("ete3 compare --taboutput -r ~/data/IES/analysis/singleGene/part.nexus.treefile -t ~/data/IES/analysis/singleGene/cluster.*.nucl.fa.renamed.treefile --unrooted > ~/data/IES/analysis/tables/singleGeneCompare.dat",0);

# infer phylogeny of concatenated genes partitioned by gene, and reusing the model test choices
#run('~/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 -spp  ~/data/IES/analysis/singleGene/part.nexus -bb 5000 > ~/data/IES/analysis/log/concatGenes.log', 1);
run('~/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 -spp  ~/data/IES/analysis/singleGene/part.nexus -redo > ~/data/IES/analysis/log/concatGenes.log', 1);


# for each gene family calculate the likelihood of the species tree and compare to the likelihood of the gene tree. Look at distribution
# TODO:
# 2. do not overright previous output add prefix(?)
# 3. create a species tree without T. thermophila for the alignments that are lacking one
# ~/tools/iqtree-1.4.2-Linux/bin/iqtree -te ~/data/IES/analysis/singleGeneOld/part.nexus.treefile -s ~/data/IES/analysis/singleGeneOld/cluster.9857.nucl.fa.renamed -redo
