#!/usr/bin/perl
use warnings;
use strict;

#prepare and run biphy

my $treePath = '/home/dsellis/data/IES_data/msas/phyldog1/results/';
my $matricesPath = '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';
my $matricesASRP = '/home/dsellis/data/IES_data/msas/asr/matrices/';
my $treeFile ='/home/dsellis/data/IES_data/msas/asr/allTrees.tre';
my $biphyPath = '/home/dsellis/tools/biphy-master/';

# read all output trees from Phyldog run
# concatenate them to one file if they have a character matrix in phylip format

opendir DH, $treePath or die "$treePath $!";
my @files = grep{/^.*.\.ReconciledTree$/} readdir(DH);
close DH;
open TR, '>'.$treeFile or die $!;
my $debug = 0;
foreach my $fileName (sort @files){
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    my $matrixF = $matricesPath.'cluster.'.$cluster.'.phy';
    if (-e $matrixF){
	open IN, $treePath.$fileName or die $!;
	my @tree = <IN>;
	close IN;
	print TR @tree;
	print $cluster,"\n";
        #copy matrix to asr/matrices
	system("cp $matrixF $matricesASRP");
    }else{
	print $cluster," has no matrix\n";
    }
    if ($debug>10){last;}
    $debug++;

}
close TR;
my $cmdl;
#$cmdl = $biphyPath."multibiphy -d $matricesASRP -t $treeFile -a ies";
$cmdl = $biphyPath."multibiphy -d $matricesASRP -t $treeFile -a -u 1 ies";
print $cmdl,"\n";
#system "$cmdl";
