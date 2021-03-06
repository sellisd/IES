#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile catdir);
use lib'.';
use functions;

# find which ies columns are within Gblocks

# read gblocks.dat

my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};
my $tablesP    = catdir($basePath, 'analysis', 'tables');

my $gblocksF = catfile($tablesP, 'gblocks.dat');
my $homcolF  = catfile($tablesP, 'homColumns.be');
my $homcolBF = catfile($tablesP, 'homColinB.be');

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
open OUT, '>', $homcolBF or die $!;
while (my $line = <HC>){
    chomp $line;
    (my $geneFamily, my $hcStart, my $hcEnd, my $ies) = split " ", $line;
    my @ies = split ",", $ies;
    if(!defined $gblocks{$geneFamily}){
	print $geneFamily, " has no Gblocks\n";
    }
    #loop through homologous IES positions
    for(my $i = 0; $i<= $#{$gblocks{$geneFamily}{'begin'}}; $i++){
	# check if within gblock
	my $bStart = $gblocks{$geneFamily}{'begin'}[$i];
	my $bEnd = $gblocks{$geneFamily}{'end'}[$i];
	if(($hcStart >= $bStart) and ($hcStart < $bEnd)
	   and(($hcEnd > $bStart) and ($hcEnd <= $bEnd))){
	    # fully within
# id geneFamily msaColumn
	    print OUT $line, "\n";
	}
    }

}
close HC;
close OUT;

