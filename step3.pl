#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;
my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.csv');
my $rbResultsP = catfile($homeD, 'data/IES/analysis/asr/');

# directories with rb output
my $rbRun1 = catfile($rbResultsP, 'run1');
my $rbRun2 = catfile($rbResultsP, 'run2');
my $rbNodeIndexesP = catfile($rbResultsP, 'rbNodeIndexes');

# Rev scripts
my $rbRun1Rev = catfile($rbRun1, 'asr1.Rev');
my $rbRun2Rev = catfile($rbRun2, 'asr2.Rev');

my $nr = getNotation($notationF);
my $asrRevF = '/home/dsellis/projects/IES/src/asr.Rev';


my $cmdl = 'Rscript --vanilla preRev.R > /home/dsellis/data/IES/analysis/preRev.log';
run($cmdl, 1);

make_path($rbNodeIndexesP) unless -e $rbNodeIndexesP;
make_path($rbRun1) unless -e $rbRun1;
make_path($rbRun2) unless -e $rbRun2;

setOutAsr($asrRevF, $rbRun1, $rbRun1Rev, 1); # first time through print node index
setOutAsr($asrRevF, $rbRun2, $rbRun2Rev, 0);
my $rbP = '/home/dsellis/tools/revbayes/projects/cmake/rb';
run("$rbP $rbRun1Rev", 1);
run("$rbP $rbRun2Rev", 1);

run('./parseRevBayesAsr.pl -gf /home/dsellis/data/IES/analysis/asr/geneFamilies.dat -burnin 1000 -asr /home/dsellis/data/IES/analysis/asr/ -output ~/data/IES/analysis/tables/avNodeProb.dat', 1);

# make a dictionary of node Ids across software
run('./nhxNodes.pl ~/data/IES/analysis/phyldog/results/*.ReconciledTree', 1);
run("./nhxNodes.pl -rb ".$rbNodeIndexesP."/nodeIndex.*.tre", 1);

# make dictionary of node indexes linking revBayes, PHYLDOG and ape(R) notation
run("Rscript --vanilla ./nodeDictionary.R", 1);

# find age of individual IES (MRCA of group of homologous IES)
# 1. spEvents.py creates a table with nodes of trees that correspond to speciation events
run("./spEvents.py > ~/data/IES/analysis/tables/spEvents.dat", 1);
# 2. then find speciation tree node is the most recent common ancestor
run("./firstIES.py > ~/data/IES/analysis/tables/firstIES.dat", 1);
# 3. then find the most recent common ancestor
run("./iesAge.py > ~/data/IES/analysis/tables/iesAge.dat", 0);

# 4. create homIESdb with age and other information
run("./addAge.pl > ~/data/IES/analysis/iesdb/homIESdb.tab", 0);

# find per branch gain and loss events
# calculate total length (nt) of conserved blocks in alignments for each gene family

# calculate all pairs of speciation nodes on species tree
run("./nodePairs.py > ~/data/IES/analysis/tables/spNodePairs.dat", 1);
# calculate  node paths for each pair of speciation nodes
run("./nodePaths.py > ~/data/IES/analysis/tables/nodePaths.dat", 1);

# for all paths connecting speciation nodes (Nanc-N1-N2-Noffspring)
# calculate the difference in probability at each step
# sum all the positive differences and all the negative differences
run("./gainLoss.pl > ~/data/IES/analysis/tables/gainLoss.dat", 1);

# for each path calculate probability of gain and loss    


sub setOutAsr{
    # modify the basic Rev script for multiple runs, and optionally keep specific lines
    my $file = shift @_;
    my $outpath = shift @_;
    my $outfile = shift @_;
    my $optline = shift @_; # should uncomment optional lines?
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
