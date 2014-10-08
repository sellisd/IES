#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Transcript;
#make one big file for each species in genbank format with one entry per contig
#TODO : find what to do with IES kai coordinates (me i xoris IES?) na tropopoiiso to pos handle species gia na einai eukolo gia ola
# make post-processing script that merges CDSes with join

# read assembled contigs file and push sequences into a hash
# read gene, cds and protein files and push into hashes
# read gff with Tools::GFF and for each entry find the sequence
#   populate it with features and the appropriate DNA and protein sequence from the appropriate files


#species name
my $speciesAbr = 'PBI';
my $species = 'Paramecium biaurelia';
my $taxonId = 65126;
#file and paths for input
my $cds = 'working/cds.fa';
my $protein = 'working/protein.fa';
my $gene = 'working/gene.fa';
my $gff3 = 'working/pb.gff3';
my $scaffoldsF = 'working/scaffolds.fa';

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
my %transcriptFeatureH;

#open sequence files and fill hashes
my $scaffoldIn = Bio::SeqIO->new('-file' => $scaffoldsF,
				 '-format' => 'fasta');
while(my $scaffoldSeq = $scaffoldIn->next_seq){
  my $scId = $scaffoldSeq->display_id;
  $scId =~ /(scaffold_\d+)_with/ or die $scId;
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

#read gff3 file and build genbank entries
while(my $feature = $gff3In->next_feature()){ # one line at a time
  my $scaffold = $feature->seq_id(); #print out the scaffold
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
  
  if ($feature->primary_tag() eq 'gene'){
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
    my $newFeature = new Bio::SeqFeature::Generic(-start => $feature->start(),
						  -end => $feature->end,
						  -strand      => $feature->strand(),
						  -primary_tag => $feature->primary_tag().' '.($feature->get_tag_values('Parent'))[0],
						 );


    # my $newFeature = new Bio::SeqFeature::Generic(-start => $feature->start(),
    #  						  -end => $feature->end(),
    #  						  -strand => $feature->strand(),
    #  						  -primary_tag => $feature->primary_tag(),
    # 						  -tag => { gene => $feature->get_tag_values('Parent'),			   #				                           translation => $proteinH{$speciesAbr.'GNP'.$number}->seq(),
    # 							    codon_start => $feature->frame()}
    # 						 );



    #print Dumper $newFeature;
    $entriesH{$scaffold}->add_SeqFeature($newFeature);#$geneFeatureH{$feature->get_tag_values('Parent')};
  }
  
  #before going to the next gene add the features to the sequence and the write the sequence
#Cant find how to print the join opperator, I will do it in two passes. Print the file and then pass second time to make exons into joins
  # my $newFeature = new Bio::SeqFeature::Gene::Transcript();
  # my $exon1 = new Bio::SeqFeature::Gene::Exon(-start => 1,
  # 						 -end => 3);
  # my $exon2 = new Bio::SeqFeature::Gene::Exon(-start => 5,
  # 						 -end => 8);
  # $newFeature->add_exon($exon1);
  # $newFeature->add_exon($exon2);
  # my $seqO = Bio::Seq->new('-display_id' => "test",
  # 			 '-format' =>'genbank');
  # $seqO->add_SeqFeature($newFeature);
  # $data_out->write_seq($seqO);
  # print Dumper $newFeature;
  # die;
}

foreach my $entry (sort keys %entriesH){
  $data_out->write_seq($entriesH{$entry});
}

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
