#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
# create fasta file with IES sequence

#read ies gff
my $input = '/home/dsellis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51_with_ies.gff3';
my $output = '>/home/dsellis/data/IES_data/ptetraurelia/ies.fa';

my $outS = Bio::SeqIO->new('-file' => $output,
			   '-format' => 'fasta');
my $inS = Bio::Tools::GFF->new('-file' => $input,
			       '-gff_version'=> 3);
while(my $feature = $inS->next_feature()){
    my $id = ($feature->get_tag_values('ID'))[0];
    my $seq =  ($feature->get_tag_values('sequence'))[0];
    my $seqO = Bio::Seq->new('-seq'      => $seq,
			     '-id'       => $id,
			     '-alphabet' => 'dna');
    $outS->write_seq($seqO);			 
}
