#!usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

# rename T. thermophila proteins
my $in = '/home/dsellis/data/IES/thermophila/gene/T_thermophila_June2014.protein.fa';

my $in  = Bio::SeqIO->new(-file => $inputFile,
			  -format => 'Fasta');
my $out = Bio::AlignIO->new(-file => '>'.$outputFile,
			    -format => 'Fasta');
while ( my $seq = $in->next_seq() ) {

    $out->write_seq($seq);
}
 
