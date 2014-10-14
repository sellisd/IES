#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Transcript;
use Bio::Location::Split;

#everything should be in MAC coordinates
# IES_junction 122..123
# idx_ref = NAME (unique)
# note "in CDS
# note stopCodon\
#  score
# logika ola einai se MAC coordinates
#make genbank file for all the rest first
#then for ies and then find in which genes they are in
#make one big file for each species in genbank format with one entry per contig

#species name
my $speciesAbr = 'PSEX';
my $species;
my $taxonId;
#file and paths for input
my $cds;
my $protein;
my $gene;
my $gff3;
my $scaffoldsF;
my $outputFile;
my $iesgffF;
#load defaults
if ($speciesAbr eq 'PBI'){
    $species = 'Paramecium biaurelia';
    $taxonId = 65126;
    $cds = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.cds.fa';
    $protein = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.protein.fa';
    $gene = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.gene.fa';
    $gff3 = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.gff3';
    $scaffoldsF = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_assembly_v1.fasta';
    $outputFile = '/Users/diamantis/data/IES_data/pbiaurelia/PBI.gnbk';
    $iesgffF = '/Users/diamantis/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3';
}elsif($speciesAbr eq 'PSEX'){
    $species = 'Paramecium sexaurelia';
    $taxonId = 65128;
    $cds = '/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.cds.fa';
    $protein = '/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.protein.fa';
    $gene = '/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.gene.fa';
    $gff3 = '/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3';
    $scaffoldsF = '/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta';
    $outputFile = '/Users/diamantis/data/IES_data/psexaurelia/PSEX.gnbk';
    $iesgffF = '/Users/diamantis/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3';
}elsif($speciesAbr eq 'PTET'){
    $species = 'Paramecium tetraurelia';
    $taxonId = 5888;
    $cds = '/Users/diamantis/data/IES_data/ptetraurelia_mac_51_annotation_v2.0.4.cds.fa';
    $protein = '/Users/diamantis/data/IES_data/ptetraurelia_mac_51_annotation_v2.0.4.protein.fa';
    $gene = '/Users/diamantis/data/IES_data/ptetraurelia_mac_51_annotation_v2.0.4.gene.fa';
    $gff3 = '/Users/diamantis/data/IES_data/ptetraurelia_mac_51_annotation_v2.0.4.gff3';
    $scaffoldsF = '/Users/diamantis/data/IES_data/ptetraurelia_mac_51.fa';
    $outputFile = '/Users/diamantis/data/IES_data/ptetraurelia/PTET.gnbk';
    $iesgffF = '/Users/diamantis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3';
}else{
    die "unknown species";
}

#output file
my $data_out = Bio::SeqIO->new('-file' => '>'.$outputFile,
				 '-format' => 'genbank');

#define hashes to be filled with sequences
my %cdsH;
my %proteinH;
my %geneH;
my %scaffoldH;

#make hash for output
my %entriesH;
my %geneFeatureH;
my %CDSH;

#open sequence files and fill hashes
my $scaffoldIn = Bio::SeqIO->new('-file' => $scaffoldsF,
				 '-format' => 'fasta');
while(my $scaffoldSeq = $scaffoldIn->next_seq){
  my $scId = $scaffoldSeq->display_id;
  $scId =~ /(scaffold_\d+)/ or die $scId; #Everything is in MAC coordinates (without IES)
  $scaffoldH{$1} = $scaffoldSeq;
}

my $cdsIn = Bio::SeqIO->new('-file' => $cds,
 			    '-format' => 'fasta');
while(my $cdsSeq = $cdsIn->next_seq){
  $cdsH{$cdsSeq->display_id} = $cdsSeq;
}

my $proteinIn = Bio::SeqIO->new('-file' => $protein,
 			    '-format' => 'fasta');
while(my $proteinSeq = $proteinIn->next_seq){
  $proteinH{$proteinSeq->display_id} = $proteinSeq;
}

my $geneIn = Bio::SeqIO->new('-file' => $gene,
 			    '-format' => 'fasta');
while(my $geneSeq = $geneIn->next_seq){
  $geneH{$geneSeq->display_id} = $geneSeq;
}

my $gff3In = Bio::Tools::GFF->new('-file' => $gff3,
				  '-gff_version' => 3);

my $curGene;
my $scaffold;
my $prevScaffold;
#read gff3 file and build genbank entries
while(my $feature = $gff3In->next_feature()){ # one line at a time
  $prevScaffold = $scaffold;
  $scaffold = $feature->seq_id(); #print out the scaffold
  if(defined($entriesH{$scaffold})){

  }else{
    #prepare sequence for this entry
    $entriesH{$scaffold} = Bio::Seq->new('-display_id' => $scaffold,
					 '-format' =>'genbank');
    my $sequence = $scaffoldH{$scaffold}->seq();
    $entriesH{$scaffold}->desc($scaffold);#definition
    $entriesH{$scaffold}->alphabet('dna');
    $entriesH{$scaffold}->seq($sequence);
    my $length = $scaffoldH{$scaffold}->length();
    my $sourceFeat = new Bio::SeqFeature::Generic(
						  -primary_tag => 'source',
						  -start       => 1,
						  -end         => $length,
						  -tag         => {organism => $species,
								   db_xref  => 'taxon:'.$taxonId}
						  #organism taxonID number
						 );
  $entriesH{$scaffold}->add_SeqFeature($sourceFeat);
  }
  my @id = $feature->get_tag_values('ID');
  my $number = deparseNumber($id[0]);

  if ($feature->primary_tag() eq 'gene'){
    if(defined($curGene)){
      if(defined($CDSH{$curGene})){
	$entriesH{$prevScaffold}->add_SeqFeature($CDSH{$curGene});
	$curGene = $number; #set the current gene
      }
    }else{
      $curGene = $number;
    }
    #build features
    #find all entries that are in the same scaffold
    # print $geneH{'PBIGNG'.$1}; #the new gene
    #print $feature->get_tag_values('ID'); die;
    $geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
							 -end         => $feature->end(),
							 -strand      => $feature->strand(),
							 -primary_tag => $feature->primary_tag(),
							 -tag => {gene     => $id[0]});
    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
  }elsif($feature->primary_tag() eq 'CDS'){
    my $parent = ($feature->get_tag_values('Parent'))[0];
    $parent = deparseNumber($parent);
    if(defined($CDSH{$parent})){
      my $location = $CDSH{$parent}->location(); #if already seen find the location of CDS
      # and add the new one to the list
      $location -> add_sub_Location(Bio::Location::Simple->new(-start  => $feature->start(),
							       -end    => $feature->end(),
							       -strand => $feature->strand()));
      $CDSH{$parent}->set_attributes(-location => $location); #update location
    }else{
      my $location = Bio::Location::Split->new();
      $location -> add_sub_Location(Bio::Location::Simple->new(-start  => $feature->start(),
							       -end    => $feature->end(),
							       -strand => $feature->strand()));
      
      $CDSH{$parent} = new Bio::SeqFeature::Generic(-location => $location,
						   -strand => $feature->strand(),
						   -primary_tag => $feature->primary_tag(),
						    -tag => { gene => $feature->get_tag_values('Parent'),
							      translation => $proteinH{$speciesAbr.'GNP'.$number}->seq()}
	  );
    }
  }elsif($feature->primary_tag() eq 'exon'){
      $geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
							   -end         => $feature->end(),
							   -strand      => $feature->strand(),
							   -primary_tag => $feature->primary_tag(),
							   -tag => {gene     => ($feature->get_tag_values('Parent'))[0],
								    codon_start => $feature->frame()});
      $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
  }
}

$entriesH{$prevScaffold}->add_SeqFeature($CDSH{$curGene}); #add the last CDS

foreach my $entry (sort keys %entriesH){
    $data_out->write_seq($entriesH{$entry});
}

$gff3In->close;

#call postProcess.pl for adding the IES information
print "adding IES information postProcess.pl:\n";
print "  ./postProcess.pl $iesgffF $outputFile\n";
exec "./postProcess.pl $iesgffF $outputFile";

sub deparseNumber{
    my $string = shift @_;
    $string =~ /.{6}(\d+)/; #extract the number regardless of whether it is a gene, exon etc
    return $1;
}
