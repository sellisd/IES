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
#TODO : find what to do with IES kai coordinates (me i xoris IES?) na tropopoiiso to pos handle species gia na einai eukolo gia ola

# read assembled contigs file and push sequences into a hash
# read gene, cds and protein files and push into hashes
# read gff with Tools::GFF and for each entry find the sequence
#   populate it with features and the appropriate DNA and protein sequence from the appropriate files

#species name
my $speciesAbr = 'PBI';
my $species = 'Paramecium biaurelia';
my $taxonId = 65126;
#file and paths for input
my $cds = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.cds.fa';
my $protein = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.protein.fa';
my $gene = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.gene.fa';
my $gff3 = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.gff3';
my $scaffoldsF = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_assembly_v1.fasta';

#output file
my $data_out = Bio::SeqIO->new('-file' => '>output.gnbk',
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
  
  $id[0] =~ /.{6}(\d+)/; #extract the number regardless of whether it is a gene, exon etc
  my $number = $1;
  # if(defined($curGene)){
  #   if($curGene ne $feature->get_tag_values('Parent')){ #if a new gene starts add the previous CDS
  #  #   $entriesH{$scaffold}->add_SeqFeature($CDSH{$id[0]});
  #   }
  #}
  if ($feature->primary_tag() eq 'gene'){
    if(defined($curGene)){
      # print $curGene," "; 
      # print keys %CDSH,"\n";
      # print $CDSH{$curGene},"\n";
      # die;
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
    #print as: CDS gene name    start..end
    # my $newFeature = new Bio::SeqFeature::Generic(-start => $feature->start(),
    # 						  -end => $feature->end,
    # 						  -strand      => $feature->strand(),
    # 						  -primary_tag => $feature->primary_tag().' '.($feature->get_tag_values('Parent'))[0],
    # 						 );
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
						   -tag => { gene => $feature->get_tag_values('Parent'),			   				                           translation => $proteinH{$speciesAbr.'GNP'.$number}->seq()}
						  );
    }
  }elsif($feature->primary_tag() eq 'exon'){
    $geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
							 -end         => $feature->end(),
							 -strand      => $feature->strand(),
							 -primary_tag => $feature->primary_tag(),
							 -tag => {gene     => ($feature->get_tag_values('Parent'))[0],
								  codon_start => $feature->frame()});
    #TODO add also reading frame or starting codon
    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
  }
}

$entriesH{$prevScaffold}->add_SeqFeature($CDSH{$curGene}); #add the last CDS
 
foreach my $entry (sort keys %entriesH){
  $data_out->write_seq($entriesH{$entry});
}

#print Dumper %CDSH;
$gff3In->close;
exit;
#read the various formats and combine information


# my $cdsIn = Bio::SeqIO->new('-file' => $cds,
# 			    '-format' => 'fasta');
# my $proteinIn = Bio::SeqIO->new('-file' => $protein,
# 			    '-format' => 'fasta');
# my $geneIn = Bio::SeqIO->new('-file' => $gene,
# 			    '-format' => 'fasta');
# my $gff3In = Bio::Tools::GFF->new('-file' => $gff3,
# 				  '-gff_version' => 3);

# #read first the gene

# #! the sequencies might not alwyas be in the same order, make sure to check in the end
# my $data_out = Bio::SeqIO->new('-file' => '>output.gnbk',
# 			       '-format' => 'genbank');

# while (my $seq = $gff3In->next_seq){
#   $data_out->write_seq($seq);
# }
# exit;
# #write feature with CDS and annotation

# while (my $geneseq = $geneIn->next_seq){
#   my $cdsseq = $cdsIn->next_seq;
#   my $proteinseq = $proteinIn->next_seq;
# #make feature object
#   # my $feat = new Bio::SeqFeature::Generic('-start' => $cdsIn,
#   # 					  '-stop' => ,
#   # 					  '-strand' => ,
#   # 					  '-primary_tag' =);

#   $data_out->write_seq($cdsseq);
# #die;
# }

# exit;

sub deparseNumber{
  my $string = shift @_;
  $string =~ /.{6}(\d+)/; #extract the number regardless of whether it is a gene, exon etc
  return $1;
}
