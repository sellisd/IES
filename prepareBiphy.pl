#!/usr/bin/perl
use warnings;
use strict;

#prepare and run biphy

my $treePath = '/home/dsellis/data/IES_data/msas/phyldog1/results/';
my $matricesPath = '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';
my $treeFile ='/home/dsellis/data/IES_data/msas/alignments/allTrees.tre';
my $biphyPath = '/home/dsellis/tools/biphy-master/';

# read all output trees from Phyldog run
# concatenate them to one file if they have a character matrix in phylip format

opendir DH, $treePath or die "$treePath $!";
my @files = grep{/.*\.ReconciledTree$/} readdir(DH);
close DH;
open TR, '>'.$treeFile or die $!;
foreach my $fileName (@files){
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    my $matrixF = $matricesPath.'cluster.'.$cluster.'.phy';
    if (-e $matrixF){
	open IN, $treePath.$fileName or die $!;
	my @tree = <IN>;
	close IN;
	print TR @tree;
	print $cluster,"\n";
    }else{
	print $cluster," has no matrix\n";
    }
}
close TR;
my $cmdl = $biphyPath."multibiphy -d $matricesPath -t $treeFile -a ies";
print $cmdl,"\n";
system "$cmdl";
