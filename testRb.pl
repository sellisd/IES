#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use File::Spec::Functions qw(catfile);

# fast test to determine which gene families have too low starting probability
my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};

my $rbP      = '/home/dsellis/tools/revbayes-1.0.0/projects/cmake/rb';
my $asrTestF = "./asrTestTemplate.Rev";
my $asr      = $ARGV[0];

open GF, catfile($basePath, 'analysis', $asr, 'geneFamilies.dat') or die $!;
chomp(my @gfs = <GF>);
close GF;

foreach my $clusterId (@gfs){
    my $outF = catfile($basePath, 'tempdat', 'testrb.'.$asr.'.'.$clusterId);
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

