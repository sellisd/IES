#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

# fast test to determine which gene families have too low starting probability

my $rbP = '/home/dsellis/tools/revbayes-1.0.0/projects/cmake/rb';
my $pm        = Parallel::ForkManager->new(7);
my $asrTestF = "/home/dsellis/projects/IES/src/asrTestTemplate.Rev";
my $asr = $ARGV[0];

open GF, '/home/dsellis/data/IES/analysis/'.$asr.'/geneFamilies.dat' or die $!;
chomp(my @gfs = <GF>);
close GF;

foreach my $clusterId (@gfs){
    my $outF = '/home/dsellis/data/IES/tempdat/testrb.'.$asr.'.'.$clusterId;
    open OUT, '>', $outF or die $!;
    open IN, $asrTestF or die $!;
    while(my $line = <IN>){
	$line =~ s/clusterS = 1000/clusterS = $clusterId/;
	$line =~ s/asr1/$asr/;
	print OUT $line;
    }
    close OUT;
    close IN;
    print "$rbP $outF\n";
}

