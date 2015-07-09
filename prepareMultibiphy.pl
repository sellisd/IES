#!/usr/bin/perl
use warnings;
use strict;
use File::Path qw(make_path);
#use Parallel::ForkManager;

#my $pm = Parallel::ForkManager->new(7);

#prepare and run biphy

my $treePath = '/home/dsellis/data/IES_data/msas/phyldog/results/';
my $matricesPath = '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';
my $matricesASRP = '/home/dsellis/data/IES_data/msas/asr/matrices/';
my $treeFile ='/home/dsellis/data/IES_data/msas/asr/allTrees.tre';
my $biphyPath = '/home/dsellis/tools/biphy-master/';

# read all output trees from Phyldog run
# concatenate them to one file if they have a character matrix in phylip format

make_path($matricesASRP) unless -d $matricesASRP;

open TR, '>'.$treeFile or die $!;

opendir DH, $treePath or die "$treePath $!";
my @files = grep{/^.*.\.ReconciledTree$/} readdir(DH);
close DH;
#my $debug = 1;

foreach my $fileName (sort @files){
#    my $pid = $pm->start and next;
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    
#    $matricesASRP = '/home/dsellis/data/IES_data/msas/asr'.$cluster.'/matrices/';
#    $treeFile ='/home/dsellis/data/IES_data/msas/asr'.$cluster.'/allTrees.tre';
 
    my $matrixF = $matricesPath.'cluster.'.$cluster.'.phy';
    if (-e $matrixF){
	open IN, $treePath.$fileName or die $!;
	my @tree = <IN>;
	close IN;
	print TR @tree;
	print $cluster,"\n";
        #copy matrix to asr/matrices
	# open MI, $matrixF or die $!;
	# open MO,'>'.$matricesASRP or die $!;
	# while(my $line = <MI>){
	#     print MO $line;
	# }
	# close MI;
	# close MO;
	system("cp $matrixF $matricesASRP");
    }else{
	print $cluster," has no matrix\n";
#	$pm->finish;
    }

#    while(! -s "ies.$cluster.trace"){
#	#wait
#    }
#    my $string = `ps -aux|grep ies.$cluster`;
#    (my $user, my $PID) = (split " ", $string)[0,1];
#    kill 'KILL', $PID;
#    $pm->finish;
#    $debug++;
#    die if $debug>10;
}
close TR;
my $cmdl;
# #$cmdl = $biphyPath."multibiphy -d $matricesASRP -t $treeFile -a ies";
$cmdl = $biphyPath."multibiphy -d $matricesASRP -t $treeFile -a -u 1 ies";
print $cmdl,"\n";
system "$cmdl";
