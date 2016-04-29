#!/usr/bin/perl
use warnings;
use strict;

# find which ies columns are within Gblocks
# read gblocks.dat

my $gblocksF = '/home/dsellis/data/IES/analysis/msas/filtered/gblocks.dat';
my $homcolF = '/home/dsellis/data/IES/analysis/tables/homColumns.be';

my %gblocks;
open GB, $gblocksF or die $!;
while (my $line = <GB>){
    chomp $line;
    (my $geneFamily, my $blockStart, my $blockEnd) = split " ", $line;
    $blockStart--; # from 1 to 0 indexed!
    if(defined($gblocks{$geneFamily})){
	push @{$gblocks{$geneFamily}{'begin'}}, $blockStart;
	push @{$gblocks{$geneFamily}{'end'}}, $blockEnd;
    }else{
	$gblocks{$geneFamily} = {'begin' => [$blockStart],
				 'end'   => [$blockEnd]
	}
    }
}
close GB;

# read homColumns.be
open HC, $homcolF or die $!;
while (my $line = <HC>){
    chomp $line;
    (my $geneFamily, my $hcStart, my $hcEnd, my $ies) = split " ", $line;
    my @ies = split ",", $ies;
    if(!defined $gblocks{$geneFamily}){
	die $geneFamily;
    }
    #loop through homologous IES positions
    for(my $i = 0; $i<= $#{$gblocks{$geneFamily}{'begin'}}; $i++){
	# check if within gblock
	my $bStart = $gblocks{$geneFamily}{'begin'}[$i];
	my $bEnd = $gblocks{$geneFamily}{'end'}[$i];
	if(($hcStart >= $bStart) and ($hcStart < $bEnd)
	   and(($hcEnd > $bStart) and ($hcEnd <= $bEnd))){
	    # fully within
	    print $line, "\n";
	}
    }

}
close HC;

