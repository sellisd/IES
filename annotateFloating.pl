#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

#read tab file (output from iesInfo.pl) with ies information and check if an IES is floating

my $tabF = $ARGV[0];
my $floatF = $tabF;
my $bedF = $tabF;
$floatF =~ s/\.tab/.float/ or die $!;
$bedF =~ s/\.tab/.be/ or die $!;
open IN, $tabF or die $!;
open FL, '>'.$floatF or die $!;
open BD, '>'.$bedF or die $!;
my $header = readline(IN);
my @outheader = ('id', 'scaffold', 'altSeqNo', 'start', 'end', 'isFloating', 'startLocs', 'endLocs', 'noAltLocs');
print FL join("\t", @outheader),"\n";
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
    push @startLocs, $start;
    push @endLocs, $end;
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
    }
    my @rowEntry = ($id, $scaffold, $altSeqNo, $start, $end, $isFloating, join(',', @startLocs), join(',', @endLocs), $#startLocs);
    print FL join("\t", @rowEntry),"\n";

    # make also a .bed file
    for(my $i = 0; $i <= $#startLocs; $i++){
	my @bedRowEntry = ($scaffold, $startLocs[$i], $endLocs[$i], $id, $isFloating);
	print BD join("\t", @bedRowEntry), "\n";
    }
}
close IN;
close FL;
close BD;
