#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
my $dataPath = '/Users/diamantis/data/IES_data/ptetraurelia/';
my $inputF = 'PTET.IES.gnbk';
my $gnbkIn = Bio::SeqIO->new('-file' => $dataPath.$inputF,
			     '-format' => 'genbank');
my $cdsout = $inputF;
$cdsout =~ s/.IES.gnbk/.CDS.bed/;
my $iesout = $inputF;
$iesout =~ s/.IES.gnbk/.IES.bed/;

open CDSOUT, '>'.$dataPath.$cdsout or die $!;
open IESOUT, '>'.$dataPath.$iesout or die $!;

while(my $seqO = $gnbkIn->next_seq()){
    foreach my $featureO ($seqO->get_SeqFeatures()){
#	print $featureO->primary_tag(),"\n";die;
	if($featureO->primary_tag() eq 'CDS'){
	    my @geneName = $featureO->get_tag_values('gene');
	    my $locationO = $featureO->location();
 	    foreach my $locations ($locationO->each_Location()){
 		my $start = $locations->start;
 		my $end = $locations->end;
		print CDSOUT "$geneName[0] $start $end\n";
	    }
	}
	if($featureO->primary_tag() eq 'IES_junction'){
	    my $locationO = $featureO->location();
	    my @id = $featureO->get_tag_values('id');
	    my $start = $locationO->start;
	    my $end = $locationO->end;
	    print IESOUT "$start $end\n";
	}
    }
}
		
# make bed files for CDS and IES from genbank file
close CDSOUT;
close IESOUT;
