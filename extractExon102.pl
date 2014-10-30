#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
my @species = qw/Pbi Pte Pse/;
my $PF = Bio::SeqIO->new('-file' => '>/Users/diamantis/data/IES_data/working/exon102P.fa',
			   '-format' => 'fasta');
my $NF = Bio::SeqIO->new('-file' => '>/Users/diamantis/data/IES_data/working/exon102N.fa',
			   '-format' => 'fasta');
my $tabOutF = '>/Users/diamantis/data/IES_data/working/exon102.tab';

#find the 102 exons and print a summary file with
#id start stop length coding frame
#and two fasta files one with nucleotide and the aminoacid sequence
open TO, $tabOutF or die $!;
foreach my $species (@species){
    my $gene;
    if($species eq 'Pbi'){
	$gene = '/Users/diamantis/data/IES_data/pbiaurelia/Pbi.gnbk';
    }elsif($species eq 'Pte'){
	$gene = '/Users/diamantis/data/IES_data/ptetraurelia/Pte.gnbk';
    }elsif($species eq 'Pse'){
	$gene = '/Users/diamantis/data/IES_data/psexaurelia/Pse.gnbk';
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
		my $codon_start = ($featureO->get_tag_values('codon_start'))[0];
		my $id = 'P_'.$species->species().'_'.$seqO->accession_number().'_'.$gene.'_'.$start.'_'.$end;
		if($length == 102){
		    my $newNSeq = Bio::Seq->new('-seq' => $seqO->subseq($start,$end),
						'-id' => $id,
						'-alphabet' => 'dna'
			);
		    my $translation = $newSeq->translate(-frame=>$codon_start,
							 -codontable_id = > 6);
		    my $newPSeq = Bio::Seq->new('-seq' => $translation,
						'-id' => $id,
						'-alphabet' => 'protein'
			);
		    $PF->write_seq($newPSeq);
		    $NF->write_seq($newNSeq);
		    print TO "$id\t$start\t$stop\t$length\t$codon_start\n";
		}
	    }
	}
    }
}
close TO;
