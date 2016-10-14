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
my $rbP = '/home/dsellis/tools/revbayes/projects/cmake/rb'; # revBayes path

my $nr = getNotation($notationF);
my $asrRevF = '/home/dsellis/projects/IES/src/asr.Rev';

# prepare revBayes runs
my $cmdl = 'Rscript --vanilla preRev.R > /home/dsellis/data/IES/analysis/preRev.log';
run($cmdl, 0);

for(my $asrRun = 1; $asrRun<=3; $asrRun++){
    my $rbResultsP = catfile($homeD, 'data/IES/analysis/asr'.$asrRun.'/');
    # directories with rb output
    my $rbRun1 = catfile($rbResultsP, 'run1');
    my $rbRun2 = catfile($rbResultsP, 'run2');
    my $rbNodeIndexesP = catfile($rbResultsP, 'rbNodeIndexes');
    # Rev scripts
    my $rbRun1Rev = catfile($rbRun1, 'asr1.Rev');
    my $rbRun2Rev = catfile($rbRun2, 'asr2.Rev');
    make_path($rbResultsP) unless -e $rbResultsP;
    make_path($rbNodeIndexesP) unless -e $rbNodeIndexesP;
    make_path($rbRun1) unless -e $rbRun1;
    make_path($rbRun2) unless -e $rbRun2;

    setOutAsr($asrRevF, $rbRun1, $rbRun1Rev, 1, $rbNodeIndexesP); # first time through print node index
    setOutAsr($asrRevF, $rbRun2, $rbRun2Rev, 0, $rbNodeIndexesP);
    run("$rbP $rbRun1Rev", 0);
    run("$rbP $rbRun2Rev", 0);

    run('./parseRevBayesAsr.pl -gf /home/dsellis/data/IES/analysis/asr'.$asrRun.'/geneFamilies.dat -burnin 1000 -asr /home/dsellis/data/IES/analysis/asr'.$asrRun.'/ -output ~/data/IES/analysis/tables/avNodeProb'.$asrRun.'.dat', 0);

    # make a dictionary of node Ids across software
    run('./nhxNodes.pl ~/data/IES/analysis/phyldogT'.$asrRun.'/results/*.ReconciledTree', 0);
    run("./nhxNodes.pl -rb ".$rbNodeIndexesP."/nodeIndex.*.tre", 0);

}
# make dictionary of node indexes linking revBayes, PHYLDOG and ape(R) notation
run("Rscript --vanilla ./nodeDictionary.R", 0);

# find age of individual IES (MRCA of group of homologous IES)
# 1. spEvents.py creates a table with nodes of trees that correspond to speciation events
run("./spEvents.py -p /home/dsellis/data/IES/analysis/phyldogT1/results/ -g /home/dsellis/data/IES/analysis/asr1/geneFamilies.dat > ~/data/IES/analysis/tables/spEvents1.dat", 0);
run("./spEvents.py -p /home/dsellis/data/IES/analysis/phyldogT2/results/ -g /home/dsellis/data/IES/analysis/asr2/geneFamilies.dat > ~/data/IES/analysis/tables/spEvents2.dat", 0);
run("./spEvents.py -p /home/dsellis/data/IES/analysis/phyldogT3/results/ -g /home/dsellis/data/IES/analysis/asr3/geneFamilies.dat > ~/data/IES/analysis/tables/spEvents3.dat", 0);

# 2. then find speciation tree node is the most recent common ancestor
run("./firstIES.py 1 > ~/data/IES/analysis/tables/iesAge1.dat", 0);
run("./firstIES.py 2 > ~/data/IES/analysis/tables/iesAge2.dat", 0);
run("./firstIES.py 3 > ~/data/IES/analysis/tables/iesAge3.dat", 0);

# 3. create homIESdb with age and other information
run("./addAge.pl 1 > ~/data/IES/analysis/iesdb/homIESdb1.tab", 0);
run("./addAge.pl 2 > ~/data/IES/analysis/iesdb/homIESdb2.tab", 0);
run("./addAge.pl 3 > ~/data/IES/analysis/iesdb/homIESdb3.tab", 0);

# find per branch gain and loss events
# calculate total length (nt) of conserved blocks in alignments for each gene family

# calculate all pairs of speciation nodes on species tree
run("./nodePairs.py /home/dsellis/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree > ~/data/IES/analysis/tables/spNodePairs1.dat", 0);
run("./nodePairs.py /home/dsellis/data/IES/analysis/phyldogT2/results/OutputSpeciesTree_ConsensusNumbered.tree > ~/data/IES/analysis/tables/spNodePairs2.dat", 0);
run("./nodePairs.py /home/dsellis/data/IES/analysis/phyldogT3/results/OutputSpeciesTree_ConsensusNumbered.tree > ~/data/IES/analysis/tables/spNodePairs3.dat", 0);

# calculate  node paths for each pair of speciation nodes
run("./nodePaths.py 1 > ~/data/IES/analysis/tables/nodePaths1.dat", 0);
run("./nodePaths.py 2 > ~/data/IES/analysis/tables/nodePaths2.dat", 0);
run("./nodePaths.py 3 > ~/data/IES/analysis/tables/nodePaths3.dat", 0);

# for all paths connecting speciation nodes (Nanc-N1-N2-Noffspring)
# calculate the difference in probability at each step
# sum all the positive differences and all the negative differences
run("./gainLoss.pl 1 > ~/data/IES/analysis/tables/gainLoss1.dat", 0);
run("./gainLoss.pl 2 > ~/data/IES/analysis/tables/gainLoss2.dat", 0);
run("./gainLoss.pl 3 > ~/data/IES/analysis/tables/gainLoss3.dat", 0);

sub setOutAsr{
    # modify the basic Rev script for multiple runs, and optionally keep specific lines
    my $file = shift @_;
    my $outpath = shift @_;
    my $outfile = shift @_;
    my $optline = shift @_; # should uncomment optional lines?
    my $rbNodeIndexesP = shift @_;
    open IN, $file or die $!;
    open OUT, '>', $outfile or die $!;
    while(my $line = <IN>){
	$line =~ s/OUTPATH/$outpath/g;
	if($optline){
	    $line =~ s/^#OPTLINE(.*)$/$1/;
	    $line =~ s/rbNodeIndexesP/$rbNodeIndexesP/;
	}
	print OUT $line;
    }
    close OUT;
    close IN;
}
