#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

# read fasta file with scaffolds

my $file = $ARGV[0];
# for each sequence
my $dataIn = Bio::SeqIO->new('-file'   => $file,
			     '-format' => 'fasta');
print "scaffold\tlength\tgc\n";
while(my $scaffold = $dataIn->next_seq){
    my $id = $scaffold->id();
    my $length = $scaffold->length();
    my $seq = $scaffold->seq();
    my $gc = $seq =~ tr/[gcGC]//;
    print $id, "\t", $length, "\t", $gc/$length;
    print "\n";
}


