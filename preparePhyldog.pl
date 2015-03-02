#!/usr/bin/perl
use warnings;
use strict;
use File::Path qw(make_path);

my $inputFasta = '/home/dsellis/data/IES_data/msas/alignments/aln/';
my $inputCharM = '/home/dsellis/data/IES_data/msas/alignments/charMat/'; #use only alignment files for which we have a character matrix
my $inputTreesF = '/home/dsellis/data/IES_data/msas/tree/trees.tre';


my $outputPath = '/home/dsellis/data/IES_data/msas/phyldog/';
my $fastaPath = 'FastaFiles/';
my $linkPath = 'LinkFiles/';
my $treePath = 'TreeFiles/';

my $speciesTree = '(Paramecium_caudatum:0.3,(Paramecium_sexaurelia:0.15,(Paramecium_tetraurelia:0.1,Paramecium_biaurelia:0.1)):0.2):0.01;';

#make required folders
make_path($outputPath.$fastaPath) unless -d $outputPath.$fastaPath;
make_path($outputPath.$linkPath) unless -d $outputPath.$linkPath;
make_path($outputPath.$treePath) unless -d $outputPath.$treePath;

open ST, '>'.$outputPath.'speciesTree.tre' or die;
print ST $speciesTree."\n";
close ST;

#read all trees in memory
my %treesH;
open TRE, $inputTreesF or die $!;
while(my $line = <TRE>){
    chomp $line;
    (my $treeName, my $treeString) = split " ", $line;
    my $cluster = substr($treeName,8,);
    $treesH{$cluster} = $treeString;
}
close TRE;
#find which clusters have character matrices
opendir(DH,$inputCharM) or die $!;
my @charMats = grep {/cluster\.\d+\.dat/} readdir(DH);
foreach my $file (@charMats){
    $file =~ /cluster\.(\d+)\.dat/;
    my $cluster = $1;
    print $cluster,"\n";
    my $fastaFileSource = $inputFasta.'cluster.'.$cluster.'.nucl.fa';
    my $fastaFileTarget = $outputPath.$fastaPath.'cluster.'.$cluster.'.fasta';
    my $linkFile = $outputPath.$linkPath.'cluster.'.$cluster.'.link';
    my $treeFile = $outputPath.$treePath.'cluster.'.$cluster.'.tre';
    open FAI, $fastaFileSource or die $!;
    open FAO, '>'.$fastaFileTarget or die $!;
    my %linkH;
    while(my $line = <FAI>){
	print FAO $line;
	chomp $line;
	if (substr($line,0,1) eq '>'){
	    my $geneName = substr($line,1);
	    my $speciesName = &gene2species($geneName);
	    if(defined($linkH{$speciesName})){
		push @{$linkH{$speciesName}}, $geneName;
	    }else{
		$linkH{$speciesName} = [$geneName];
	    }
	}
    }
    close FAI;
    close FAO;
    open LIN, '>'.$linkFile or die $!;
    foreach my $speciesName (sort keys %linkH){
	print LIN $speciesName,':';
	print LIN join(';', @{$linkH{$speciesName}}),"\n";	
    }
    close LIN;
    open TROUT, '>'.$treeFile or die $!; # <')(((((<
    if(defined($treesH{$cluster})){
	print TROUT $treesH{$cluster};
    }else{
# 3 gene tree!
	}
    close TROUT;
}
close DH;

sub gene2species{
#get species name from gene name
    my $string = shift @_;
    my $speciesName;
    my $abr = substr($string,0,4);
    if($abr eq 'PCAU'){
	$speciesName = 'Paramecium_caudatum';
    }elsif($abr eq 'PSEX'){
	$speciesName = 'Paramecium_sexaurelia';
    }elsif($abr eq 'PTET'){
	$speciesName = 'Paramecium_tetraurelia';
    }elsif($abr eq 'PBIA'){
	$speciesName = 'Paramecium_biaurelia';
    }else{
	die "unknown name $abr";
    }
}
#prepare appropriate directories and files for a phyldog run

#make 3 directories
# FastaFiles
#   cluster.X.fasta
# LinkFiles
#   cluster.X.link
#      speciesName:geneName1;geneName2;
# TreeFiles
#
