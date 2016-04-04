#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Transcript;
use Bio::Location::Split;
use Getopt::Long;
use Bio::Species;
use File::Spec::Functions qw(catfile);
use File::HomeDir;
use lib'.';
use functions;

my $help;
my $species3abr;
my $dataPath;
my $floating = 1; #by default filter floating IES
my $usage = <<HERE;

make genebank file from gff3 files
    usage makeGeneBank.pl -species Pab -datapath PATH -floating 1

consistent species abreviations, these are not the same with the abbreviations provided in the sequence files
 tth Tetrahymena thermophila
 pca Paramecium caudatum
 ppr Paramecium primaurelia
 pbi Paramecium biaurelia
 pte Paramecium tetraurelia
 ppe Paramecium pentaurelia
 pse Paramecium sexaurelia
 poc Paramecium octaurelia
 ptr Paramecium tredecaurelia
 pso Paramecium sonneborni

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'datapath=s' => \$dataPath,
			      'species=s' => \$species3abr));
die $usage if $help;
die $usage unless defined($species3abr);

my $homeD = File::HomeDir->my_home;
$dataPath = catfile($homeD, 'data/IES/'); #default for local run
#$dataPath = '/pandata/IES_data'; #default for cluster

#everything should be in MAC coordinates
# IES_junction 122..123
# idx_ref = NAME (unique)
# note "in CDS
# note stopCodon\
#  score

#make genbank file for all the rest first
#then for ies and then find in which genes they are in
#make one big file for each species in genbank format with one entry per contig
my @lineage;
if( $species3abr eq 'ppr' or
    $species3abr eq 'pbi' or
    $species3abr eq 'pte' or
    $species3abr eq 'ppe' or
    $species3abr eq 'pse' or
    $species3abr eq 'poc' or
    $species3abr eq 'ptr' or
    $species3abr eq 'pso' or
    $species3abr eq 'pca'
    ){
    @lineage  = ('Eukaryota','Alveolata','Ciliophora','Intramacronucleata','Oligohymenophorea','Peniculida','Parameciidae','Paramecium');
}elsif($species3abr eq 'tth'){
    @lineage  = ('Eukaryota','Alveolata','Ciliophora','Intramacronucleata', 'Oligohymenophorea', 'Hymenostomatida','Tetrahymenina', 'Tetrahymenidae', 'Tetrahymena');
}else{
    die "unknown species";
}

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
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.tab');

open N, $notationF or die $!;
my $header = readline(N);
my %notation;
while(my $line = <N>){
    chomp $line;
    my $ar = split "\t", $line;
    (my $abbreviation, my $datapath, my $binomial, my $taxId, my $geneGff, my $cdsF, my $protF, my $geneF, my $MacF, my $iesGff, my $annotation, my $prefix) = split "\t", $line;
    $notation{$abbreviation} = {
	'datapath'      => $datapath,
	'binomial'	=> $binomial,
	'taxId'	        => $taxId,
	'geneGff'	=> $geneGff,
	'cdsF'	        => $cdsF,
	'protF'	        => $protF,
	'geneF'	        => $geneF,
	'MacF'	        => $MacF,
	'iesGff'	=> $iesGff,
	'annotation'    => $annotation,
	'prefix'        => $prefix
    };
}
close N;

#initialize prefix names:
my $prefixR = initF();

unless (defined($notation{$species3abr})){
    die "unknown species abbreviation $species3abr";
}

$species    = $notation{$species3abr}{'binomial'};
$taxonId    = $notation{$species3abr}{'taxId'};
$speciesAbr = $notation{$species3abr}{'abbreviation'};
$dataPath   = $notation{$species3abr}{'datapath'};
$cds        = catfile($dataPath, $notation{$species3abr}{'cdsF'});
$protein    = catfile($dataPath, $notation{$species3abr}{'protF'});
$gene       = catfile($dataPath, $notation{$species3abr}{'geneF'});
$gff3       = catfile($dataPath, $notation{$species3abr}{'geneGff'});
$scaffoldsF = catfile($dataPath, $notation{$species3abr}{'MacF'});
$outputFile = catfile($dataPath, 'analysis/gnbk/'.$species3abr.'.gnbk');
$iesgffF    = catfile($dataPath, $notation{$species3abr}{'iesGff'});

push @lineage, $species;
@lineage = reverse(@lineage);
my $speciesO = Bio::Species->new(-classification => \@lineage);

#output file
my $data_out = Bio::SeqIO->new('-file' => '>'.$outputFile,
			       '-format' => 'genbank');

#define hashes to be filled with sequences
my %cdsH;
my %proteinH;
my %geneH;
my %scaffoldH;
my %geneIdNameH;

#make hash for output
my %entriesH;
my %geneFeatureH;
my %CDSH;
my %otherH; # hash for all the rest of the features

my %rna2geneH; # hash linking mrna to gene ids

my %scaffoldElementsH; # for each scaffold list with all elements

#open sequence files and fill hashes
print "read scaffolds from $scaffoldsF\n";
my $scaffoldIn = Bio::SeqIO->new('-file' => $scaffoldsF,
				 '-format' => 'fasta');
while(my $scaffoldSeq = $scaffoldIn->next_seq){
    my $scId = $scaffoldSeq->display_id;
    $scaffoldH{$scId} = $scaffoldSeq;
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
    if(defined($proteinH{$id})){
	die;
    }else{
	$proteinH{$id} = $proteinSeq;
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

while(my $feature = $gff3In->next_feature()){ # one line at a time
    my $scaffold = $feature->seq_id(); # scaffold
    if(defined($entriesH{$scaffold})){
	
    }else{
	#prepare sequence for this entry
	$entriesH{$scaffold} = Bio::Seq->new('-display_id' => $scaffold,
					     '-format' =>'genbank',
					     '-accession_number' => $scaffold);
	$entriesH{$scaffold}->species($speciesO);
	my $sequence = $scaffoldH{$scaffold}->seq();
	$entriesH{$scaffold}->desc($scaffold); # definition
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
    if(defined($scaffoldElementsH{$scaffold})){
	push @{$scaffoldElementsH{$scaffold}},$id[0];
    }else{
	$scaffoldElementsH{$scaffold} = [$id[0]];
    }
    if ($feature->primary_tag() eq 'gene'){
	$geneIdNameH{$id[0]}=($feature->get_tag_values('Name'))[0];
	#find all entries that are in the same scaffold
	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
							     -end         => $feature->end(),
							     -strand      => $feature->strand(),
							     -primary_tag => $feature->primary_tag(),
							     -tag => {gene     => $id[0]});
    }elsif($feature->primary_tag() eq 'CDS'){
	my $parent = ($feature->get_tag_values('Parent'))[0];
	my $geneId;# = ($feature->get_tag_values('Parent'))[0];
	if(defined($rna2geneH{$parent})){
	    $geneId = $rna2geneH{$parent};
	}else{
	    die;
	}
	if(defined($CDSH{$geneId})){
	    my $location = $CDSH{$geneId}->location(); #if already seen find the location of CDS
	    # and add the new one to the list
	    $location -> add_sub_Location(Bio::Location::Simple->new(-start  => $feature->start(),
								     -end    => $feature->end(),
								     -strand => $feature->strand()));
	    $CDSH{$geneId}->set_attributes(-location => $location); #update location
	}else{
	    my $location = Bio::Location::Split->new();
	    $location -> add_sub_Location(Bio::Location::Simple->new(-start  => $feature->start(),
								     -end    => $feature->end(),
								     -strand => $feature->strand()));
	    my $proteinId = &getProtName($geneId, $species3abr);
	    $CDSH{$geneId} = new Bio::SeqFeature::Generic(-location => $location,
							  -strand => $feature->strand(),
							  -primary_tag => $feature->primary_tag(),
							  -tag => { gene => $geneId,
								    codon_start => $feature->frame(),
								    translation => $proteinH{$proteinId}->seq()
							  }
		);
	}
    }elsif(($feature->primary_tag() eq 'mRNA')
	   or ($feature->primary_tag() eq 'ncRNA')){ #link exons to genes through mRNA
	my $id = ($feature->get_tag_values('ID'))[0];
	my $parent = ($feature->get_tag_values('Parent'))[0];
	$rna2geneH{$id} = $parent;
	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start        => $feature->start(),
							     -end          => $feature->end(),
							     -strand       => $feature->strand(),
							     -primary_tag  => $feature->primary_tag(),
							     -tag => {gene => $parent,
							     }
	    );
    }elsif($feature->primary_tag() eq 'exon'){
	my $parent = ($feature->get_tag_values('Parent'))[0];
	if(defined($rna2geneH{$parent})){
	    $parent = $rna2geneH{$parent};
	}else{
	    die "$parent";
	}
	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start        => $feature->start(),
							     -end          => $feature->end(),
							     -strand       => $feature->strand(),
							     -primary_tag  => $feature->primary_tag(),
							     -tag => {gene => $parent,
							     }
	    );
#	    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
    }else{
#get all tags, if it has parent and the parent is a gene
	my $found = 0;
	my @tags = $feature->get_all_tags();
	foreach my $tag(@tags){
	    if ($tag eq 'Parent'){
		$found = 1;
	    }
	}
	if($found){
	    my $parent = ($feature->get_tag_values('Parent'))[0];
	    if(defined($rna2geneH{$parent})){
		$parent = $rna2geneH{$parent};
	    }else{
		die "$parent";
	    }
	    $geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start        => $feature->start(),
								 -end          => $feature->end(),
								 -strand       => $feature->strand(),
								 -primary_tag  => $feature->primary_tag(),
								 -tag => {gene => $parent,
								 }
		);
#		$entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
	}else{
	    $otherH{$feature->primary_tag()} = 1;	    
	}
    }
}

print "printing output\n";

#add features to scaffolds
foreach my $scaffold(sort keys %scaffoldElementsH){
    foreach my $element (@{$scaffoldElementsH{$scaffold}}){
	if(defined($geneFeatureH{$element})){
	    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$element});
	    if(defined($CDSH{$element})){
		$entriesH{$scaffold}->add_SeqFeature($CDSH{$element});
	    }
	}
    }
    $data_out->write_seq($entriesH{$scaffold});
}

print '  did not include in features: ',join(',', keys(%otherH)),"\n";

$gff3In->close;

#call postProcess.pl for adding the IES information
print "adding IES information postProcess.pl:\n";
print "  ./postProcess.pl -species $species3abr $iesgffF $outputFile\n";
exec "./postProcess.pl -species $species3abr $iesgffF $outputFile";

sub getProtName{
    # find name of protein from name of gene.
    # if paramecium
    #   protein name P<->G
    # if tetrahymena 
    #   protein name is gene name
    my $geneId = shift @_;
    my $proteinId;
    if($species3abr eq 'th'){
	$proteinId = $geneIdNameH{$geneId};
    }else{
	$proteinId = gene2prot($geneId, $prefixR);
    }
    if(defined($proteinId)){
	return $proteinId;
    }else{
	die "unable to find protein name for: $geneId";
    }
}
