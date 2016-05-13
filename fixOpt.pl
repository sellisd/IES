#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);

# correct option files created by prepareData.py
my $inD = $ARGV[0];
my $outD = $ARGV[1];

opendir DH, $inD or die $!;
my @files = grep {/^.*\.opt$/} readdir(DH);
foreach my $file (@files){
    open IN, catfile($inD, $file) or die $!;
    open OUT, '>', catfile($outD, $file) or die $!;
    while(my $line = <IN>){
	$line =~ s/^(input\.sequence\.file=\/pandata\/sellis\/phyldog\/aln\/)opt\/(\d+)\.opt/$1$2.fasta/;
	print OUT $line;
    }
    close OUT;
    close IN;
}
close DH;
