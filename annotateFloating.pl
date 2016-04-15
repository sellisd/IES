#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Getopt::Long;
my $help;
my $tabF;
my $bedF;
my $floatF;
my $usage = <<HERE;

read tab file (output from iesInfo.pl) with ies information and check if an IES is floating

usage annotateFloating.pl [OPTIONS]

where OPTIONS can be:

    -tabF:   input tab file
    -floatF: output file with floating annotation
    -bedF:   output file in bed format

HERE

die $usage unless(GetOptions(
		      'help|?'   => \$help,
		      'tabF=s'   => \$tabF,
		      'bedF=s'   => \$bedF,
		      'floatF=s' => \$floatF
		  ));

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
    	my $start = $startLocs[$i];
   	    $start--; # gff3 is [start, end] while in .bed is (start, end]
		my @bedRowEntry = ($scaffold, $start, $endLocs[$i], $id, $isFloating);
		print BD join("\t", @bedRowEntry), "\n";
    }
}
close IN;
close FL;
close BD;
