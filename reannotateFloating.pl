#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# check if IES require reannotation. Read output from annotate floating and for each location check if multiple IES are expected to be present

my @inFiles = @ARGV;
foreach my $anotF (@inFiles){
#    print $anotF," \n";
    my %loc; #scaffold.location => $id
    open IN, $anotF or die $!;
    while(my $line = <IN>){
	chomp $line;
	(my $id, my $scaffold, my $altSeqNo, my $start, my $end, my $isFloating, my $startLocs, my $endLocs, my $startLocNo) = split "\t", $line;
	my @startLocs = split(',', $startLocs);
	foreach my $floatLoc (@startLocs){
	    if(defined($loc{$scaffold.'.'.$floatLoc})){
		#another IES is annotated in the same region!!!
		&printab($anotF, $id, $isFloating, $loc{$scaffold.'.'.$floatLoc});#	    $loc{$scaffold.'.'.$startLocs} = $id;
	    }else{
		$loc{$scaffold.'.'.$floatLoc} = $id.' '.$isFloating;
	    }
	}
    }
    close IN;
}
