#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

#read tab file and check for floating

my $tabF = $ARGV[0];
open IN, $tabF or die $!;
my $header = readline(IN);
&printab('id', 'scaffold', 'altSeqNo', 'start', 'end', 'isFloating', 'startLocs', 'endLocs', 'noAltLocs');
while(my $line = <IN>){
    chomp $line;
    (my $id, my $scaffold, my $altSeqNo, my $start, my $end, my $upstreamFlank, my $downstreamFlank, my $length, my $front, my $back, my $sequence) = split "\t", $line;
    my $float;
    eval{
	$float = &isFloating($sequence, $upstreamFlank, $downstreamFlank);
    };
    if($@){
	print $@;
	print $line,"\n";
	die;
    }
    my $isFloating;
    my @startLocs;
    my @endLocs;
    if (ref($float) eq 'ARRAY'){
	#is floating
	$isFloating = 1;
	foreach my $displacement(@$float){
	    push @startLocs, $displacement + $start;
	    push @endLocs, $displacement + $end;
	}
    }else{
	if($float eq '0'){
	    $isFloating = 0;
	}elsif($float eq '+?'){
	    $isFloating = 'NA';
	}elsif($float eq '-?'){
	    $isFloating = 'NA';
	}
	push @startLocs, $start;
	push @endLocs, $end;
    }
    &printab($id, $scaffold, $altSeqNo, $start, $end, $isFloating, join(',', @startLocs), join(',', @endLocs), $#startLocs);
}
close IN;
