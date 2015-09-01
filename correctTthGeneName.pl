#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
my $fileGnbk = '/home/dsellis/data/IES_data/tthermophila/Tth.r.seq';
my $fileGff3 = '/home/dsellis/data/IES_data/tthermophila/T_thermophila_June2014.gff3';
my $out = '/home/dsellis/data/IES_data/tthermophila/Tth.r.seq.1';

my $data_out = Bio::SeqIO->new('-file' => '>'.$out,
			       '-format' => 'genbank');
# read gff file make hash with id and name
my $gff3In = Bio::Tools::GFF->new('-file' => $fileGff3,
				  '-gff_version' => 3);
my %idsH;
my %namesH;
while(my $feature = $gff3In->next_feature()){
    next unless($feature->primary_tag eq 'gene');
    my @tags = $feature->get_all_tags;
    my @id;
    my @name;
    my $found = 0;
    foreach my $tag (@tags){
	if($tag eq 'ID'){
	    @id = $feature->get_tag_values($tag);
	    $found = 1;
	}
	if($tag eq 'Name'){
	    @name = $feature->get_tag_values($tag);
	    $found = 1;
	}
    }
    next unless $found == 1;
    die if(defined($idsH{$id[0]})); #make sure they are unique
    die if(defined($namesH{$name[0]}));
    $idsH{$id[0]}=$name[0];
    $namesH{$name[0]} = $id[0];
}
# read genbank file and replace gene name
# use Data::Dumper;
# print Dumper %namesH;die;

my $gnbkIn = Bio::SeqIO->new('-file' => $fileGnbk,
			     '-format' => 'genbank');

while(my $seqO = $gnbkIn->next_seq()){
    my $id = $seqO->display_id;
    for my $featO ($seqO->get_SeqFeatures){
	my @tags = $featO->get_all_tags;
	foreach my $tag (@tags){
	    if($tag eq 'gene'){
		my @geneId = $featO->get_tag_values('gene');
		$featO->remove_tag('gene');
		$featO->add_tag_value('gene',$idsH{$geneId[0]});
	    }
	}
    }
    $data_out->write_seq($seqO);
}
