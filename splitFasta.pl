#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $dataPath = '/pandata/sellis/';
my $fastaOutPath = '/pandata/sellis/msas/fasta/';


#load groupings in memory
print "loading silix output in memory\n";
my %hash;
my $silixOutput = $dataPath.'working/silix.output';

open IN, $silixOutput or die $!;
while (my $line = <IN>){
    chomp $line;
    (my $group, my $id) = split " ", $line;
    if(defined($hash{$group})){
	push @{$hash{$group}}, $id;
    }else{
	$hash{$group}=[$id];
    }
}
close IN;
