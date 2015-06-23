#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

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
	print $cluster,"\n";
	my $inblocksF = 'cluster.'.$cluster.'.inblocks.bed';
	# read inblocksF make arrays
	my @blockStart;
	my @blockEnd;
	open GB, $pathFiltered.$gbF or die $!;
	while(my $line = <GB>){
	    chomp $line;
	    (my $cluster, my $start, my $end) = split " ", $line;
	    push @blockStart, $start;
	    push @blockEnd, $end;
	}
	close GB;
	my @iesStart;
	my @iesEnd;
	open BD, $pathFiltered.$bdF or die $!;
	while(my $line = <BD>){
	    chomp $line;
	    (my $cluster, my $start, my $end, my $id) = split " ", $line;
	    push @iesStart, $start;
	    push @iesEnd, $end;
	}
	close BD;
	my $inOneR = whichInOne(\@iesStart,\@iesEnd,\@blockStart,\@blockEnd);
	
	open OUT, '>'.$pathFiltered.$inblocksF or die $!;
	foreach my $in1 (sort{$a<=>$b} keys %$inOneR){
	    print OUT "$cluster\t$iesStart[$in1]\t$iesEnd[$in1]\t$in1\n";
	}
	close OUT;
    }
}

# find which IES are partially overlaping
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my $inBlocksF = 'cluster.'.$cluster.'.inblocks.bed';
    if(-s $pathFiltered.$inBlocksF){ # if there are IES in blocks
	my $ies2mergeCMD = 'bedtools merge -o collapse -c 4 -i '.$pathFiltered.$inBlocksF.' > '.$pathFiltered.'cluster.'.$cluster.'.2merge.bed';
	print "$ies2mergeCMD\n";
	system "$ies2mergeCMD";
    }
}
