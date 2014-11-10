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
use Getopt::Long;
use Bio::Species;
my $help;
my $species3abr;
my $dataPath;
my $usage = <<HERE;

make genebank file from gff3 files
usage makeGeneBank.pl -species Pab -datapath PATH

#consistent species abreviations, these are not the same with the abbreviations provided in the sequence files
# Ppr Paramecium primaurelia
# Pbi Paramecium biaurelia
# Pte Paramecium tetraurelia
# Ppe Paramecium pentaurelia
# Pse Paramecium sexaurelia

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'datapath=s' => \$dataPath,
			      'species=s' => \$species3abr));
die $usage if $help;
$dataPath = '/Users/diamantis/data/IES_data/'; #default for local run
#$dataPath = '/pandata/IES_data'; #default for cluster

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
my @lineage  = ('Eukaryota','Alveolata','Ciliophora','Intramacronucleata','Oligohymenophorea','Peniculida','Parameciidae','Paramecium');

my $speciesAbr;
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
if ($species3abr eq 'Pbi'){
    $species = 'Paramecium biaurelia';
    $taxonId = 65126;
    $speciesAbr = 'PBIA.V1_4.1.';
    $dataPath = $dataPath.'pbiaurelia/';
    $cds = $dataPath.'pbiaurelia_V1-4_annotation_v2.0.cds.fa';
    $protein = $dataPath.'pbiaurelia_V1-4_annotation_v2.0.protein.fa';
    $gene = $dataPath.'pbiaurelia_V1-4_annotation_v2.0.gene.fa';
    $gff3 = $dataPath.'pbiaurelia_V1-4_annotation_v2.0.gff3';
    $scaffoldsF = $dataPath.'biaurelia_V1-4_assembly_v1.fasta';
    $outputFile = $dataPath.'Pbi.gnbk';
    $iesgffF = $dataPath.'internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3';
}elsif($species3abr eq 'Pse'){
    $species = 'Paramecium sexaurelia';
    $taxonId = 65128;
    $speciesAbr = 'PSEX.AZ8_4.';
    $dataPath = $dataPath.'psexaurelia/';
    $cds = $dataPath.'psexaurelia_AZ8-4_annotation_v2.0.cds.fa';
    $protein = $dataPath.'psexaurelia_AZ8-4_annotation_v2.0.protein.fa';
    $gene = $dataPath.'psexaurelia_AZ8-4_annotation_v2.0.gene.fa';
    $gff3 = $dataPath.'psexaurelia_AZ8-4_annotation_v2.0.gff3';
    $scaffoldsF = $dataPath.'sexaurelia_AZ8-4_assembly_v1.fasta';
    $outputFile = $dataPath.'Pse.gnbk';
    $iesgffF = $dataPath.'internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3';
}elsif($species3abr eq 'Pte'){
    $species = 'Paramecium tetraurelia';
    $taxonId = 5888;
    $speciesAbr = 'PTET.51.1.';
    $dataPath = $dataPath.'ptetraurelia/';
    $cds = $dataPath.'ptetraurelia_mac_51_annotation_v2.0.cds.fa';
    $protein = $dataPath.'ptetraurelia_mac_51_annotation_v2.0.protein.fa';
    $gene = $dataPath.'ptetraurelia_mac_51_annotation_v2.0.gene.fa';
    $gff3 = $dataPath.'ptetraurelia_mac_51_annotation_v2.0.gff3';
    $scaffoldsF = $dataPath.'ptetraurelia_mac_51.fa';
    $outputFile = $dataPath.'Pte.gnbk';
    $iesgffF = $dataPath.'internal_eliminated_sequence_PGM_IES51.pt_51.gff3';
}else{
    die "unknown species";
}
push @lineage, $species;
@lineage = reverse(@lineage);
my $speciesO = Bio::Species->new(-classification => \@lineage);

#output file
my $data_out = Bio::SeqIO->new('-file' => '>'.$outputFile,
			       '-format' => 'genbank');
## $data_out->write_seq($sequencedbg);
## die;
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
print "read scaffolds from $scaffoldsF\n";
my $scaffoldIn = Bio::SeqIO->new('-file' => $scaffoldsF,
				 '-format' => 'fasta');
while(my $scaffoldSeq = $scaffoldIn->next_seq){
  my $scId = $scaffoldSeq->display_id;
  $scId =~ /(scaffold_?[_\d]+)/ or die $scId; #Everything is in MAC coordinates (without IES) scaffold names are not consistent across species
  $scaffoldH{$1} = $scaffoldSeq;
}

print "read CDS from $cds\n";
my $cdsIn = Bio::SeqIO->new('-file' => $cds,
 			    '-format' => 'fasta');
while(my $cdsSeq = $cdsIn->next_seq){
  $cdsH{$cdsSeq->display_id} = $cdsSeq;
}
print "read proteins from $protein\n";
my $proteinIn = Bio::SeqIO->new('-file' => $protein,
 			    '-format' => 'fasta');
while(my $proteinSeq = $proteinIn->next_seq){
    my $id = $proteinSeq->display_id;
    my $idNo = deparseNumber($id);
    if(defined($proteinH{$idNo})){
	die;
    }else{
	$proteinH{$idNo} = $proteinSeq;
    }
}
print "read genes from $gene\n";
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
print "parse gff3 features\n";
#read gff3 file and build genbank entries
while(my $feature = $gff3In->next_feature()){ # one line at a time
    $prevScaffold = $scaffold;
    $scaffold = $feature->seq_id(); #print out the scaffold
    if(defined($entriesH{$scaffold})){
	
    }else{
	#prepare sequence for this entry
	$entriesH{$scaffold} = Bio::Seq->new('-display_id' => $scaffold,
					     '-format' =>'genbank',
					     '-accession_number' => $scaffold);
	$entriesH{$scaffold}->species($speciesO);
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
    $id[0] =~ /$speciesAbr.+?(\d+):?/;
    my $number = $1;
    if ($feature->primary_tag() eq 'gene'){
	if(defined($curGene)){
	    if(defined($CDSH{$curGene})){
		$entriesH{$prevScaffold}->add_SeqFeature($CDSH{$curGene});
		$curGene = $number; #set the current gene
	    }else{
		$curGene = $number; #for a gene without CDS
	    }
	}else{
	    $curGene = $number;
	}
	#build features
	#find all entries that are in the same scaffold
	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
							     -end         => $feature->end(),
							     -strand      => $feature->strand(),
							     -primary_tag => $feature->primary_tag(),
							     -tag => {gene     => $id[0]});
	$entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
    }elsif($feature->primary_tag() eq 'CDS'){
	my $parent = ($feature->get_tag_values('Parent'))[0];
	$parent = deparseNumber($parent);
	my $geneId = ($feature->get_tag_values('Parent'))[0];
	$geneId =~s/(.*)T(\d+)/$1G$2/;
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
							  -tag => { gene => $geneId,
								    codon_start => $feature->frame(),
								    translation => $proteinH{$number}->seq()}
		);
	}
    }elsif($feature->primary_tag() eq 'exon'){
     	my $geneId = ($feature->get_tag_values('Parent'))[0];
     	$geneId =~s/(.*)T(\d+)/$1G$2/;
     	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
     							     -end         => $feature->end(),
     							     -strand      => $feature->strand(),
     							     -primary_tag => $feature->primary_tag(),
     							     -tag => {gene     => $geneId,
								      codon_start => $feature->frame()});
     	$entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
    }else{
     	my $geneId = ($feature->get_tag_values('Parent'))[0];
     	$geneId =~s/(.*)T(\d+)/$1G$2/;
     	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
     							     -end         => $feature->end(),
     							     -strand      => $feature->strand(),
     							     -primary_tag => $feature->primary_tag(),
     							     -tag => {gene     => $geneId,
							     });
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
print "  ./postProcess.pl -species $species3abr $iesgffF $outputFile\n";
exec "./postProcess.pl -species $species3abr $iesgffF $outputFile";

sub deparseNumber{
    my $string = shift @_;
    $string =~ /$speciesAbr.*\D(\d+)/ or die; #extract the number regardless of whether it is a gene, exon etc
    return $1;
}
