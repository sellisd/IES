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
my $cmdl;
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    $cmdl .= ' -cds '.catfile($pab{'datapath'}, $pab{'cdsF'});
}
$cmdl .= ' -cds /home/dsellis/data/IES/genomicData/thermophila/gene/T_thermophila_June2014_CDS.fasta'; #add Tth

run('./prot2nucl.pl '.$cmdl.' ~/data/IES/analysis/msas/filtered/*.aln.fa', 1);
run('./prot2nucl.pl -noterm'.$cmdl.' ~/data/IES/analysis/singleGene/*.aln.fa', 1);

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
	my $cmdl = "$iqtreeB -s ".catfile($singleGeneP, $file).
#	' -m TESTNEWONLY -b 100';
	    ' -st CODON6'.
	    ' -m TESTNEW -redo';
#	    ' -m TESTNEWONLY -redo';
	run($cmdl, 1);
	$pm->finish;
    }
    $pm->wait_all_children;
}

my $bmF = '/home/dsellis/data/IES/analysis/tables/bestModels.tab';
run('./bestModel.pl -nex ~/data/IES/analysis/singleGene/part.nexus -table '.$bmF.' ~/data/IES/analysis/singleGene/cluster.*.nucl.fa.renamed', 1);

run("ete3 compare --taboutput -r ~/data/IES/analysis/singleGene/part.nexus.treefile -t ~/data/IES/analysis/singleGene/cluster.*.nucl.fa.renamed.treefile --unrooted > ~/data/IES/analysis/tables/singleGeneCompare.dat",1);

# infer phylogeny of concatenated genes partitioned by gene, and reusing the model test choices
#run('~/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 -spp  ~/data/IES/analysis/singleGene/part.nexus -bb 5000 > ~/data/IES/analysis/log/concatGenes.log', 1);
#run('~/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 -spp  ~/data/IES/analysis/singleGene/part.nexus -b 100 > ~/data/IES/analysis/log/concatGenes.log -redo', 0);
run('~/tools/iqtree-omp-1.4.2-Linux/bin/iqtree-omp -nt 7 -spp  ~/data/IES/analysis/singleGene/part.nexus -redo > ~/data/IES/analysis/log/concatGenes.log', 1);


# for each gene family calculate the likelihood of the species tree and compare to the likelihood of the gene tree. Look at distribution
# TODO:
# 2. do not overright previous output add prefix(?)
# 3. create a species tree without T. thermophila for the alignments that are lacking one
# ~/tools/iqtree-1.4.2-Linux/bin/iqtree -te ~/data/IES/analysis/singleGeneOld/part.nexus.treefile -s ~/data/IES/analysis/singleGeneOld/cluster.9857.nucl.fa.renamed -redo

# find gene families with Tth
my %bmH;
open BM, $bmF or die "$! $bmF";
readline(BM); #header
while(my $line = <BM>){
    chomp $line;
    (my $file, my $geneFamily, my $hasTth, my $bestModel) = split " ", $line;
    $bmH{$geneFamily} = [$hasTth, $bestModel];
}
close BM;
#use Data::Dumper;print Dumper %bmH;die;
my @noTth;
foreach my $geneFamily (sort {$a<=$b} keys %bmH){
    my $hasTth    = ${$bmH{$geneFamily}}[0];
    next if $hasTth == 0;
    my $bestModel = ${$bmH{$geneFamily}}[1];
    my $inputFile =  '/home/dsellis/data/IES/analysis/singleGene/cluster.'.$geneFamily.'.nucl.fa.renamed';
    $cmdl = $iqtreeB.
	' -te ~/data/IES/analysis/singleGene/part.nexus.treefile'.
	' -redo'.
	' -s '.$inputFile.
	' -st CODON6'.
	' -pre  ~/data/IES/analysis/singleGene/cluster.'.$geneFamily.'.nucl.fa.renamed.lkh'.
	' -m '.$bestModel.
	' -fixbr';
    run($cmdl, 1);
    push @noTth, $inputFile;
}

$cmdl = './parseLikelihoods.pl '.join(' ', @noTth).' > ~/data/IES/analysis/tables/likelihoodComp.dat';
run($cmdl, 1);



#iq tree relevant options
 -z trees to compute the log-likelihoods (concat)
 -zb replicates for tests
 -s alignment
 -redo
 -st CODON6
 -pre prefix
 -m GY+F3X4+R2
 -fixbr
 -zb 1000 -zw -au
Analysis results written to: 
  IQ-TREE report:                /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.lkh.iqtree
  Maximum-likelihood tree:       /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.lkh.treefile
  Likelihood distances:          /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.lkh.mldist
  Evaluated user trees:          /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.lkh.trees
  Screen log file:               /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.lkh.log

create a file with both gene tree and reference (concatenated)
cat ~/data/IES/analysis/singleGene/part.nexus.treefile  /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.treefile > ~/data/IES/analysis/singleGene/example.treefile
/home/dsellis/tools/iqtree-1.4.2-Linux/bin/iqtree -z ~/data/IES/analysis/singleGene/example.treefile -redo -s /home/dsellis/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed -st CODON6 -pre  ~/data/IES/analysis/singleGene/cluster.7487.nucl.fa.renamed.lkh -m GY+F3X4+R2 -fixbr -zb 1000 -zw -au
