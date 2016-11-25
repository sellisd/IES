#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;
my $pm        = Parallel::ForkManager->new(7);
my $asrTestF = "/home/dsellis/projects/IES/src/asrTest1.Rev";
my $asr = 'asr3';

open GF, '/home/dsellis/data/IES/analysis/'.$asr.'/geneFamilies.dat' or die $!;
chomp(my @gfs = <GF>);
close GF;

foreach my $clusterId (@gfs){
    my $outF = '/home/dsellis/data/IES/tempdat/testrb'.$clusterId;
    open OUT, '>', $outF or die $!;
    open IN, $asrTestF or die $!;
    while(my $line = <IN>){
	$line =~ s/clusterS = 1000/clusterS = $clusterId/;
	$line =~ s/asr1/$asr/;
	print OUT $line;
    }
    close OUT;
    close IN;
    print "~/tools/revbayes/projects/cmake/rb $outF\n";
}

