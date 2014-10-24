#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
my @species = qw/Pbi Pte Pse/;
my $OutF = Bio::SeqIO->new('-file' => '>/Users/diamantis/data/IES_data/working/exon102.fa',
			   '-format' => 'fasta');

foreach my $species (@species){
    my $gene;
    if($species eq 'Pbi'){
	$gene = '/Users/diamantis/data/IES_data/pbiaurelia/PBI.gnbk';
    }elsif($species eq 'Pte'){
	$gene = '/Users/diamantis/data/IES_data/ptetraurelia/PTET.gnbk';
    }elsif($species eq 'Pse'){
	$gene = '/Users/diamantis/data/IES_data/psexaurelia/PSEX.gnbk';
    }else{
	die;
    }
    
    my $geneIn = Bio::SeqIO->new('-file' => $gene,
				 '-format' => 'genBank');
    while(my $seqO = $geneIn->next_seq){
	my $scaffold = $seqO->accession_number();
	my $species = $seqO->species();
	foreach my $featureO ($seqO->get_SeqFeatures()){
	    if($featureO->primary_tag() eq 'exon'){
		my $length = $featureO->length();
		my $start = $featureO->start();
		my $end = $featureO->end();
		my $gene = ($featureO->get_tag_values('gene'))[0];
		if($length == 102){
		    my $newSeq = Bio::Seq->new('-seq' => $seqO->subseq($start,$end),
					       '-id' => 'P_'.$species->species().'_'.$seqO->accession_number().'_'.$gene.'_'.$start.'_'.$end,
			);
		    $OutF->write_seq($newSeq);
		}
	    }
	}
    }
}
