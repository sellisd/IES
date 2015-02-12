#!/usr/bin/perl
use warnings;
use strict;

# run Gblocks on the nucleotide alignments

my $path = $ARGV[0];

opendir(DH, $path) or die $!;
my @files = grep { /nucl.fa/ } readdir(DH);
foreach my $file (@files){
    my $cmdl = 'Gblocks '.$path.$file.' -t="codons"';
    print $cmdl,"\n";
    system $cmdl;
}