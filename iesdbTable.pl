#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Getopt::Long;

my %merge;
my %ies;

my $help;
my $mergeF;
my $iesinF;
my $iestabF;
my $floatF;
my $usage = <<HERE;

combine the result of analyses to create the IES table for iesDB
usage iesdbTable.pl [OPTIONS]
where OPTIONS can be:
    -merge : file with IES floating on other annotated IES
    -iesin : bed file with IES overlapping CDS, intergenic regions or introns
    -iestab: file with IES sequence information annotations extracted from gff3
    -float : file with floating status annotation
    -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?'   => \$help,
			      'merge=s'  => \$mergeF,
			      'iesin=s'  => \$iesinF,
			      'float=s'  => \$floatF,
			      'iestab=s' => \$iestabF
		   ));
die $usage if $help;

# find merged ies
open M, $mergeF or die $!;
readline(M); #header
while(my $line = <M>){
    chomp $line;
    (my $abr, my $ies1, my $ies2) = split " ", $line;
    if(defined($merge{$abr})){
	push @{$merge{$abr}}, ($ies1, $ies2);
    }else{
	$merge{$abr} = [$ies1, $ies2];
    }
}
close M;

# read ies and find if floating
open FL, $floatF or die $!;
readline(FL); #header
while(my $line = <FL>){
    chomp $line;
    (my $id, my $scaffold, my $altSeqNo, my $start, my $end, my $isFloating, my $startLocs, my $endLocs, my $noAltLocs) = split " ", $line;
    $ies{$id}{'isFloating'} = $isFloating;
    $ies{$id}{'inCDS'} = 0;
    $ies{$id}{'inInter'} = 0;
    $ies{$id}{'inIntron'} = 0;
    $ies{$id}{'startLocs'} = [split(',', $startLocs)];
    $ies{$id}{'endLocs'} = [split(',', $endLocs)];
}
close FL;

# read ies that (partially or fully) overlap some element of interest
open T, $iesinF or die $!;
while(my $line = <T>){
    chomp $line;
    (my $scaffold, my $start, my $end, my $name, my $isFloating, my $inType, my $inScaffold, my $inStart, my $inEnd, my $name1, my $overlap) = split " ", $line;
    # some sanity checks
    if ($scaffold ne $inScaffold) { die $scaffold.' '.$inScaffold;}
    die if $start >= $end;
    die if $inStart > $inEnd;
    if($inType eq 'cds'){
	$ies{$name}{'inCDS'} = $name1; # set now
    }elsif($inType eq 'intergenic'){
	$ies{$name}{'inInter'}++;
    }elsif($inType eq 'intron'){
	$ies{$name}{'inIntron'}++;
    }
}
close T;

open TAB, $iestabF or die $!;
readline(TAB); # header
printab(('id', 'scaffold', 'altSeqNo', 'startLocs', 'endLocs', 'upstreamFlank', 'downstreamFlank', 'length', 'isFloating', 'merged', 'inCDS', 'inInter', 'inIntron', 'front', 'back', 'seq'));
while(my $line = <TAB>){
    chomp $line;
    (my $id, my $scaffold, my $altSeqNo, my $start, my $end, my $upstreamFlank, my $downstreamFlank, my $length, my $front, my $back, my $seq) = split " ", $line;
    my $merged = 0;
    foreach my $sp (keys %merge){
	foreach my $ies (@{$merge{$sp}}){
	    if ($ies eq $id){
		$merged = 1;
	    }
	}
    }
    my $floating = $ies{$id}{'isFloating'};
    if($floating){
	$start = join(',', @{$ies{$id}{'startLocs'}});
	$end = join(',', @{$ies{$id}{'endLocs'}});
    }
    my @ar = ($id, $scaffold, $altSeqNo, $start, $end, $upstreamFlank, $downstreamFlank, $length, $floating, $merged, $ies{$id}{'inCDS'}, $ies{$id}{'inInter'}, $ies{$id}{'inIntron'}, $front, $back, $seq);
    printab(@ar);
}
close TAB;
