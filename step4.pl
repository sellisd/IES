#!/usr/bin/perl
use File::Spec::Functions qw(catfile catdir);
use lib'.';
use functions;
use warnings;
use strict;

my $opt = loadUserOptions;
my $basePath   = $$opt{'basePath'};
my $tablesP    = catdir($basePath, 'analysis', 'tables');
my $iesdbP     = catdir($basePath, 'analysis', 'iesdb');

# check MCMC runs
run("Rscript --vanilla ./validateMCMC.R", 1);

for(my $asrRun = 1; $asrRun<=3; $asrRun++){
   my $rbResultsP     = catdir($basePath, 'analysis', 'asr'.$asrRun);
   my $rbNodeIndexesP = catdir($rbResultsP, 'rbNodeIndexes');
   my $geneFamiliesF  = catfile($rbResultsP, 'geneFamilies.dat');
   my $avNodeProbF    = catfile($tablesP, 'avNodeProb'.$asrRun.'.dat');
   run('./parseRevBayesAsr.pl'.
       ' -gf '.$geneFamiliesF.
       ' -burnin 1000'.
       ' -asr '.$rbResultsP.
       ' -output '.$avNodeProbF,
       1);
}

# make dictionary of node indexes linking revBayes, PHYLDOG and ape(R) notation
run("./nodeDict.py", 1);

# find age of individual IES (MRCA of group of homologous IES)
# 1. spEvents.py creates a table with nodes of trees that correspond to speciation events
my $p1rP           = catdir($basePath, 'analysis', 'phyldogT1', 'results');
my $p2rP           = catdir($basePath, 'analysis', 'phyldogT2', 'results');
my $p3rP           = catdir($basePath, 'analysis', 'phyldogT3', 'results');
my $geneFamilies1F = catfile($basePath, 'analysis', 'asr1', 'geneFamilies.dat');
my $geneFamilies2F = catfile($basePath, 'analysis', 'asr2', 'geneFamilies.dat');
my $geneFamilies3F = catfile($basePath, 'analysis', 'asr3', 'geneFamilies.dat');
my $spEvents1F     = catfile($tablesP, 'spEvents1.dat');
my $spEvents2F     = catfile($tablesP, 'spEvents2.dat');
my $spEvents3F     = catfile($tablesP, 'spEvents3.dat');

run('./spEvents.py -p '.$p1rP.' -g '.$geneFamilies1F.' > '.$spEvents1F, 0);
run('./spEvents.py -p '.$p2rP.' -g '.$geneFamilies2F.' > '.$spEvents2F, 0);
run('./spEvents.py -p '.$p3rP.' -g '.$geneFamilies3F.' > '.$spEvents3F, 0);

# 2. then find speciation tree node is the most recent common ancestor
my $iesAge1F = catfile($tablesP, 'iesAge1.dat');
my $iesAge2F = catfile($tablesP, 'iesAge2.dat');
my $iesAge3F = catfile($tablesP, 'iesAge3.dat');
run("./firstIES.py 1 > $iesAge1F", 0);
run("./firstIES.py 2 > $iesAge2F", 0);
run("./firstIES.py 3 > $iesAge3F", 0);

# 3. create homIESdb with age and other information
my $homIESdb1F = catfile($iesdbP, 'homIESdb1.tab');
my $homIESdb2F = catfile($iesdbP, 'homIESdb2.tab');
my $homIESdb3F = catfile($iesdbP, 'homIESdb3.tab');

run("./addAge.pl 1 > $homIESdb1F", 1);
run("./addAge.pl 2 > $homIESdb2F", 1);
run("./addAge.pl 3 > $homIESdb3F", 1);

# find per branch gain and loss events
# calculate total length (nt) of conserved blocks in alignments for each gene family

# calculate all pairs of speciation nodes on species tree
my $treeout1F = catfile($basePath, 'analysis', 'phyldogT1', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree');
my $treeout2F = catfile($basePath, 'analysis', 'phyldogT2', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree');
my $treeout3F = catfile($basePath, 'analysis', 'phyldogT3', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree');
my $spNodePairs1F = catfile($tablesP, 'spNodePairs1.dat');
my $spNodePairs2F = catfile($tablesP, 'spNodePairs2.dat');
my $spNodePairs3F = catfile($tablesP, 'spNodePairs3.dat');

run("./nodePairs.py $treeout1F > $spNodePairs1F", 1);
run("./nodePairs.py $treeout2F > $spNodePairs2F", 1);
run("./nodePairs.py $treeout3F > $spNodePairs3F", 1);

# calculate  node paths for each pair of speciation nodes
my $nodePaths1F = catfile($tablesP, 'nodePaths1.dat');
my $nodePaths2F = catfile($tablesP, 'nodePaths2.dat');
my $nodePaths3F = catfile($tablesP, 'nodePaths3.dat');
run("./nodePaths.py 1 > $nodePaths1F", 1);
run("./nodePaths.py 2 > $nodePaths2F", 1);
run("./nodePaths.py 3 > $nodePaths3F", 1);

# for all paths connecting speciation nodes (Nanc-N1-N2-Noffspring)
# calculate the difference in probability at each step
# sum all the positive differences and all the negative differences
my $gl1F = catfile($tablesP, 'gainLoss1.dat');
my $gl2F = catfile($tablesP, 'gainLoss1.dat');
my $gl3F = catfile($tablesP, 'gainLoss1.dat');
run("./gainLoss.pl 1 > $gl1F", 1);
run("./gainLoss.pl 2 > $gl2F", 1);
run("./gainLoss.pl 3 > $gl3F", 1);

my $gainLoss1F    = catfile($tablesP, 'gainLoss1.dat');
my $gainLoss2F    = catfile($tablesP, 'gainLoss2.dat');
my $gainLoss3F    = catfile($tablesP, 'gainLoss3.dat');
my $gblocksF      = catfile($tablesP, 'gblocks.dat');
my $brlenTree1F   = catfile($basePath, 'analysis', 'sgf', 'topoConstrSimple.treefile');
my $brlenTree2F   = catfile($basePath, 'analysis', 'sgf', 'concatSimple.nexus.treefile');
my $brlenTree3F   = catfile($basePath, 'analysis', 'sgf', 'concat.nexus.treefile');
my $phyldogTree1F = catfile($basePath, 'analysis', 'phyldogT1', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree');
my $phyldogTree2F = catfile($basePath, 'analysis', 'phyldogT2', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree');
my $phyldogTree3F = catfile($basePath, 'analysis', 'phyldogT3', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree');
my $output1F      = catfile($basePath, 'analysis', 'figures', 'spTree1');
my $output2F      = catfile($basePath, 'analysis', 'figures', 'spTree2');
my $output3F      = catfile($basePath, 'analysis', 'figures', 'spTree3');

run("./gainLossSum.py -r -g $gainLoss1F -b $gblocksF -l $brlenTree1F -p $phyldogTree1F -o $output1F -n Tetrahymena_thermophila", 1);
run("./gainLossSum.py -r -g $gainLoss2F -b $gblocksF -l $brlenTree2F -p $phyldogTree2F -o $output2F -n Tetrahymena_thermophila", 1);
run("./gainLossSum.py -r -g $gainLoss3F -b $gblocksF -l $brlenTree3F -p $phyldogTree3F -o $output3F -n Tetrahymena_thermophila", 1);
# summarize and prepare plots at the end of the analysis
exit(1);
#~/anaconda_ete/bin/python ./gainLossSum.py -i ~/Desktop/temp.dat -g /Volumes/WDC/data/IES/analysis/tables/gainLoss1.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/topoConstrSimple.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree1 -n Tetrahymena_thermophila
#~/anaconda_ete/bin/python ./gainLossSum.py -i ~/Desktop/temp.dat -g /Volumes/WDC/data/IES/analysis/tables/gainLoss2.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/concatSimple.nexus.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT2/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree2 -n Tetrahymena_thermophila
#~/anaconda_ete/bin/python ./gainLossSum.py -i ~/Desktop/temp.dat -g /Volumes/WDC/data/IES/analysis/tables/gainLoss3.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/concat.nexus.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT3/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree3 -n Tetrahymena_thermophila

#~/anaconda_ete/bin/python ./gainLossSum.py -i ~/Desktop/temp.dat -r -g /Volumes/WDC/data/IES/analysis/tables/gainLoss1.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/topoConstrSimple.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree1Norm -n Tetrahymena_thermophila
#~/anaconda_ete/bin/python ./gainLossSum.py -i ~/Desktop/temp.dat -r -g /Volumes/WDC/data/IES/analysis/tables/gainLoss2.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/concatSimple.nexus.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT2/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree2Norm -n Tetrahymena_thermophila
#~/anaconda_ete/bin/python ./gainLossSum.py -i ~/Desktop/temp.dat -r -g /Volumes/WDC/data/IES/analysis/tables/gainLoss3.dat -b /Volumes/WDC/data/IES/analysis/tables/gblocks.dat -l /Volumes/WDC/data/IES/analysis/sgf/concat.nexus.treefile -p /Volumes/WDC/data/IES/analysis/phyldogT3/results/OutputSpeciesTree_ConsensusNumbered.tree -o /Volumes/WDC/data/IES/analysis/figures/spTree3Norm -n Tetrahymena_thermophila

system("./pwm.py ~/data/IES/analysis/tables/consensus.dat"); # calculate position weight matrices and draw sequence logo diagrams
system("./lengthAge.pl > ~/data/IES/analysis/tables/ageLength.dat"); # combine age and length information for each IES
