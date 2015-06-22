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

my @inblocks;
#find which IES are within a Gblocks cluster
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my $bdF = 'cluster.'.$cluster.'.bed';
    my $gbF = 'cluster.'.$cluster.'.aln.fasta.gblocks';
    if(-e $pathFiltered.$gbF){ # if alignment passed the identity and seqNo filters
	my $inblocksF = 'cluster.'.$cluster.'.inblocks.bed';
	my $iesInBlocksCMD = 'bedtools intersect -wb -a '.$pathFiltered.$gbF.' -b '.$pathFiltered.$bdF.' > '.$pathFiltered.$inblocksF;
	print "$iesInBlocksCMD\n";
	system "$iesInBlocksCMD";
	push @inblocks, $inblocksF;
    }

}

# read intermediate files
# in case an IES overlaps 2 separate Gblocks remove it!
foreach my $inblocksF (@inblocks){
    my %hash; # make sure IES are located in only one Gblock
    my %lines2del;
    open IN, $pathFiltered.$inblocksF or die $!;
    my $lineCounter = 0;
    while(my $line = <IN>){
	chomp $line;
	(my $cluster1, my $start1, my $end1, my $cluster2, my $start2, my $end2, my $id) = split " ", $line;
	if(defined($hash{$id})){
	    $lines2del{$lineCounter} = 1;
	}else{
	    $hash{$id} = 1;
	}
	$lineCounter++;
    }
    close IN;
    if(%lines2del){
	rename $pathFiltered.$inblocksF, $pathFiltered.$inblocksF.'.temp' or die $!;
	open IN, $pathFiltered.$inblocksF.'.temp' or die $!;
	open OUT, '>'.$pathFiltered.$inblocksF or die $!;
	my $lineCounter = 0;
	while(my $line = <IN>){
	    if(defined($lines2del{$lineCounter})){ #skip line(s)
	    }else{
		print OUT $line;
	    }
	    $lineCounter++;
	}
	close OUT;
	close IN;
	unlink $pathFiltered.$inblocksF.'.temp' or die $!;
    }
}
die;
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
