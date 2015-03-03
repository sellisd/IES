#!/usr/bin/perl
use warnings;
use strict;

my $path = '/home/dsellis/data/IES_data/msas/alignments/aln/';
my $pathOut = '/home/dsellis/data/IES_data/msas/alignments/charMat/';
opendir(DH, $path) or die $!;
my @charMF = grep {/\.F\.dat$/} readdir(DH);

#find which IES are within a Gblocks cluster
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my $bdF = 'cluster.'.$cluster.'.bed';
    my $gbF = 'cluster.'.$cluster.'.nucl.fa.gblocks';
    my $iesInBlocksCMD = 'bedtools intersect -wb -a '.$path.$gbF.' -b '.$path.$bdF.' > '.$path.'cluster.'.$cluster.'.inblocks.bed';
    print "$iesInBlocksCMD\n";
    system "$iesInBlocksCMD";
}

# find which IES are partially overlaping
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my $inBlocksF = 'cluster.'.$cluster.'.inblocks.bed';
    if(-s $path.$inBlocksF){ # if there are IES in blocks
	my $ies2mergeCMD = 'bedtools merge -o collapse -c 7 -i '.$path.$inBlocksF.' > '.$path.'cluster.'.$cluster.'.2merge.bed';
	print "$ies2mergeCMD\n";
	system "$ies2mergeCMD";
    }
}

#or find which IES share either start or an end
