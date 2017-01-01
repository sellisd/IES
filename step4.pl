#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use lib'.';
use functions;
use warnings;
use strict;

# check MCMC runs
run("Rscript --vanilla ./validateMCMC.R", 1);
my $basePath = '/Users/dsellis/data/IES';

for(my $asrRun = 1; $asrRun<=3; $asrRun++){
   my $rbResultsP = catfile($basePath, 'analysis/asr'.$asrRun.'/');
   my $rbNodeIndexesP = catfile($rbResultsP, 'rbNodeIndexes');
   run('./parseRevBayesAsr.pl -gf ~/data/IES/analysis/asr'.$asrRun.'/geneFamilies.dat -burnin 1000 -asr ~/data/IES/analysis/asr'.$asrRun.'/ -output ~/data/IES/analysis/tables/avNodeProb'.$asrRun.'.dat', 1);
   # make a dictionary of node Ids across software
   run('./nhxNodes.pl ~/data/IES/analysis/phyldogT'.$asrRun.'/results/*.ReconciledTree', 1);
   run("./nhxNodes.pl -rb ".$rbNodeIndexesP."/nodeIndex.*.tre", 1);
}

# make dictionary of node indexes linking revBayes, PHYLDOG and ape(R) notation
run("Rscript --vanilla ./nodeDictionary.R", 1);

# find age of individual IES (MRCA of group of homologous IES)
# 1. spEvents.py creates a table with nodes of trees that correspond to speciation events
run("./spEvents.py -p ~/data/IES/analysis/phyldogT1/results/ -g ~/data/IES/analysis/asr1/geneFamilies.dat > ~/data/IES/analysis/tables/spEvents1.dat", 1);
run("./spEvents.py -p ~/data/IES/analysis/phyldogT2/results/ -g ~/data/IES/analysis/asr2/geneFamilies.dat > ~/data/IES/analysis/tables/spEvents2.dat", 1);
run("./spEvents.py -p ~/data/IES/analysis/phyldogT3/results/ -g ~/data/IES/analysis/asr3/geneFamilies.dat > ~/data/IES/analysis/tables/spEvents3.dat", 1);

# 2. then find speciation tree node is the most recent common ancestor
run("./firstIES.py 1 > ~/data/IES/analysis/tables/iesAge1.dat", 1);
run("./firstIES.py 2 > ~/data/IES/analysis/tables/iesAge2.dat", 1);
run("./firstIES.py 3 > ~/data/IES/analysis/tables/iesAge3.dat", 1);

# 3. create homIESdb with age and other information
run("./addAge.pl 1 > ~/data/IES/analysis/iesdb/homIESdb1.tab", 0);
run("./addAge.pl 2 > ~/data/IES/analysis/iesdb/homIESdb2.tab", 0);
run("./addAge.pl 3 > ~/data/IES/analysis/iesdb/homIESdb3.tab", 0);

# find per branch gain and loss events
# calculate total length (nt) of conserved blocks in alignments for each gene family

# calculate all pairs of speciation nodes on species tree
run("./nodePairs.py ~/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree > ~/data/IES/analysis/tables/spNodePairs1.dat", 0);
run("./nodePairs.py ~/data/IES/analysis/phyldogT2/results/OutputSpeciesTree_ConsensusNumbered.tree > ~/data/IES/analysis/tables/spNodePairs2.dat", 0);
run("./nodePairs.py ~/data/IES/analysis/phyldogT3/results/OutputSpeciesTree_ConsensusNumbered.tree > ~/data/IES/analysis/tables/spNodePairs3.dat", 0);

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


exit(1);
my $dbg = "
./gainLossSum.py -g ~/data/IES/analysis/tables/gainLoss1.dat -b ~/data/IES/analysis/tables/gblocks.dat -l ~/data/IES/analysis/sgf/topoConstrSimple.treefile -p ~/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree -o ~/data/IES/analysis/figures/spTree1 -n Tetrahymena_thermophila

./gainLossSum.py -g ~/data/IES/analysis/tables/gainLoss2.dat -b ~/data/IES/analysis/tables/gblocks.dat -l ~/data/IES/analysis/sgf/concatSimple.nexus.treefile -p ~/data/IES/analysis/phyldogT2/results/OutputSpeciesTree_ConsensusNumbered.tree -o ~/data/IES/analysis/figures/spTree2 -n Tetrahymena_thermophila

./gainLossSum.py -g ~/data/IES/analysis/tables/gainLoss3.dat -b ~/data/IES/analysis/tables/gblocks.dat -l ~/data/IES/analysis/sgf/concat.nexus.treefile -p ~/data/IES/analysis/phyldogT3/results/OutputSpeciesTree_ConsensusNumbered.tree -o ~/data/IES/analysis/figures/spTree3 -n Tetrahymena_thermophila

./gainLossSum.py -g /Volumes/WDC/data/IES/analysis/tables/gainLoss1.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/topoConstrSimple.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree1 -n Tetrahymena_thermophila
~/anaconda_ete/bin/python ./gainLossSum.py -g /Volumes/WDC/data/IES/analysis/tables/gainLoss2.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/concatSimple.nexus.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT2/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree2 -n Tetrahymena_thermophila
~/anaconda_ete/bin/python ./gainLossSum.py -g /Volumes/WDC/data/IES/analysis/tables/gainLoss3.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/concat.nexus.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT3/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree3 -n Tetrahymena_thermophila

/Volumes/WDC/
";
# summarize and prepare plots at the end of the analysis

system("./pwm.py ~/data/IES/analysis/tables/consensus.dat"); # calculate position weight matrices and draw sequence logo diagrams
system("./lengthAge.pl > ~/data/IES/analysis/tables/ageLength.dat"); # combine age and length information for each IES


