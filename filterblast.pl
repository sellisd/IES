#!/usr/bin/perl
use warnings;
use strict;

my $file = $ARGV[0];

my $lengthCutoff = 20;
my $identityCutoff = 95;

# blast output
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

open IN, $file or die $!;
while(my $line = <IN>){
	my @ar = split "", $line;
	next if $ar[0] eq $ar[1]; # skip self hits
	next if $ar[2] <= $identityCutoff;
	next if $ar[3] <= $lengthCutoff;
	print $line;
}
close IN;