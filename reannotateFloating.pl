#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# check if IES require reannotation. Read output from annotate floating and for each location check if multiple IES are expected to be present

my @inFiles = @ARGV;
&printab('species', 'ies1', 'ies2');
my %tomerge;
foreach my $anotF (@inFiles){
#    print $anotF," \n";
    $anotF =~ /^.*(p..)\.ies\.float$/;
    my $species = $1;
    my %loc; #scaffold.location => $id
    open IN, $anotF or die $!;
    while(my $line = <IN>){
	chomp $line;
	(my $id, my $scaffold, my $altSeqNo, my $start, my $end, my $isFloating, my $startLocs, my $endLocs, my $startLocNo) = split "\t", $line;
	my @startLocs = split(',', $startLocs);
	foreach my $floatLoc (@startLocs){
	    if(defined($loc{$scaffold.'.'.$floatLoc})){
		#another IES is annotated in the same region!!!
		&printab($species, $id, $loc{$scaffold.'.'.$floatLoc});#	    $loc{$scaffold.'.'.$startLocs} = $id;
		$tomerge{$id} = 1;
		$tomerge{$loc{$scaffold.'.'.$floatLoc}} = 1;
	    }else{
		$loc{$scaffold.'.'.$floatLoc} = $id;
	    }
	}
    }
    close IN;


}
