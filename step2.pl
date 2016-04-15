#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;

my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.tab');

open N, $notationF or die $!;
my $header = readline(N);
my %notation;
while(my $line = <N>){
    chomp $line;
(my $abbreviation, my $datapath, my $binomial, my $taxId, my $geneGff, my $cdsF, my $protF, my $geneF, my $MacF, my $iesGff, my $annotation, my $prefix) = split "\t", $line;
    $notation{$binomial} = {
	'annotation' => $annotation,
	'prefix'     => $prefix,
	'abr'        => $abbreviation
    };
}
close N;

if(1){
    system "./msaLocal.pl > msaLocal.log";
}

if(1){
    system "Rscript --vanilla ./filterProtAlign.R"
}
