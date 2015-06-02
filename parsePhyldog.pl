#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;
# print the MPR trees
# parse species tree and keep names of nodes
#parse phyldog output
# MRCA
# 1. P2      P. biaurelia
# 2. P4      P. tetraurelia
# 3. P6      P. sexaurelia
# 4. Pc      P. caudatum
# 5. T       T. thermophila
# 6. P24     P. biaurelia & P. tetraurelia
# 7. P246    P. biaurelia & P. tetraurelia & P. sexaurelia
# 8. P246c   P. biaurelia & P. tetraurelia & P. sexaurelia & P. caudatum
# 9. P246cT  P. biaurelia & P. tetraurelia & P. sexaurelia & P. caudatum & T. thermophila

#The node numbers are hardcoded (data from OutputSpeciesTree_ConsensusNumbered.tree)

my $MRCA1 = 8;  #  P. biaurelia
my $MRCA2 = 7;  #  P. tetraurelia
my $MRCA3 = 5;  #  P. sexaurelia
my $MRCA6 = 6;  #  P. biaurelia & P. tetraurelia 
my $MRCA7 = 4;  #  P. biaurelia & P. tetraurelia & P. sexaurelia


my $file = '/home/dsellis/data/IES_data/msas/phyldog/result/10691.ReconciledTree';

my $asrF = '/home/dsellis/data/IES_data/msas/alignments/asr/10691.dat';

my %hash;
open IN, $asrF or die $!;
while (my $line = <IN>){
    chomp $line;
    next if $line eq 'lower upper';
    (my $node, my $lower, my $upper) = split " ", $line;
    $hash{$node} = [$lower,$upper];
}


my $input = new Bio::TreeIO(-file => $file, -format => "nhx"); my
    $tree = $input->next_tree; for my $node ($tree->get_nodes){
	my @Ev = @{$node->{'_tags'}->{'Ev'}};
	my @S = @{$node->{'_tags'}->{'S'}};
	my @ND = @{$node->{'_tags'}->{'ND'}};
#	print $node->id,"\t";
	print "Ev:@Ev S:@S ND:@ND\n";

#    print "   ",$node->bootstrap,"\n";
}
use Data::Dumper;
###print Dumper $tree;
