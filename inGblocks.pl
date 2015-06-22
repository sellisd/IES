#!/usr/bin/perl
use warnings;
use strict;

# find which IES are within Gblocks only for alignments that have passed the (%) identity and sequence number filters
# and merges partially overlapping IES (floating)

#my $path = '/home/dsellis/data/IES_data/msas/alignments/aln/';
my $pathFiltered = '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $pathOut = '/home/dsellis/data/IES_data/msas/alignments/charMat/';

opendir(DH, $pathFiltered) or die $!;
my @charMF = grep {/\.F\.dat$/} readdir(DH);

#find which IES are within a Gblocks cluster
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my $bdF = 'cluster.'.$cluster.'.bed';
    my $gbF = 'cluster.'.$cluster.'.aln.fasta.gblocks';
    if(-e $pathFiltered.$gbF){ # if alignment passed the identity and seqNo filters
	my $iesInBlocksCMD = 'bedtools intersect -wb -a '.$pathFiltered.$gbF.' -b '.$pathFiltered.$bdF.' > '.$pathFiltered.'cluster.'.$cluster.'.inblocks.bed';
	print "$iesInBlocksCMD\n";
	system "$iesInBlocksCMD";
    }
}

# find which IES are partially overlaping
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my $inBlocksF = 'cluster.'.$cluster.'.inblocks.bed';
    if(-s $pathFiltered.$inBlocksF){ # if there are IES in blocks
	my $ies2mergeCMD = 'bedtools merge -o collapse -c 7 -i '.$pathFiltered.$inBlocksF.' > '.$pathFiltered.'cluster.'.$cluster.'.2merge.bed';
	print "$ies2mergeCMD\n";
	system "$ies2mergeCMD";
    }
}
