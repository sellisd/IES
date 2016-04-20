#!/usr/bin/perl
use warnings;
use strict;

# run Gblocks on the protein alignments

my $path = $ARGV[0];

opendir(DH, $path) or die $!;
my @files = grep { /cluster\..*\.aln\.fa$/ } readdir(DH);
foreach my $file (@files){
    my $cmdl = 'Gblocks '.$path.$file.' -t="protein"';
    print $cmdl,"\n";
    system $cmdl;
}
