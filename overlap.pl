#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Bio::Tools::GFF;
#alternative approach:
#./makeBed.pl creates bed files for IES and CDS from the IES.gnbk file
#using bedtools calculate overlap, which genes have an IES in them
# bedtools intersect -a PTET.CDS.bed -b PTET.IES.bed -wo > Pte.overlap
my %overlapH;
my $overlapF = 'Pte.overlap';
open OR, $overlapF or die $!;
while(my $line = <OR>){
    chomp $line;
    my @ar = split " ", $line;
    my $scaffold = $ar[0];
    my $CDSiD = $ar[3];
    my $IESiD = $ar[7];
    if(defined($overlapH{$scaffold}{$CDSiD})){
	push @{$overlapH{$scaffold}{$CDSiD}},$IESiD;
    }else{
	$overlapH{$scaffold}{$CDSiD} = [$IESiD];
    }
}
close OR;
print Dumper %overlapH; die;

# parse overlap and load CDSs overlapping an IES
#$hash{$scaffold}{$CDS.id => [IES.id,...]}

#find where in proteins IES are located
#read IES locations
#read CDS file
#calculate overlap
#print protein index and IES name
# read genbank file
#foreach scaffold
#  foreach CDS loop through segments

#read IES locations make $hash{scaffold}->{start => id}
my %iesH;
my $iesgffF = '/Users/diamantis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3';

my $iesgffO = Bio::Tools::GFF->new('-file' => $iesgffF,
				  '-gff_version' => 3);

while(my $feature = $iesgffO->next_feature()){
    my $scaffold = $feature->seq_id(); #print out the scaffold
    my $start = $feature->start();
    my $name = $feature->get_tag_values('ID');
    $iesH{$scaffold}{$start}=$name;
}

my $gnbkIn = Bio::SeqIO->new('-file' => $ARGV[0],
			     '-format' => 'genbank');
while(my $seqO = $gnbkIn->next_seq()){
    foreach my $featureO ($seqO->get_SeqFeatures()){
	if($featureO->primary_tag() eq 'CDS'){
	    my $locationO = $featureO->location();
	    my $index = 0;
	    foreach my $locations ($locationO->each_Location()){
		my $start = $locations->start;
		my $end = $locations->end;
#if there is an IES between start end

		#print gene IES name index+ (IES location - CDS start+1)
		$index += ($end-$start+1) #length of each CDS
	    }
	}
    }
}
