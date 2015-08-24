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

my $help;
my $species3abr;
my $dataPath;
my $floating = 1; #by default filter floating IES
my $usage = <<HERE;

make genebank file from gff3 files
    usage makeGeneBank.pl -species Pab -datapath PATH -floating 1

#consistent species abreviations, these are not the same with the abbreviations provided in the sequence files
# Tth Tetrahymena thermophila
# Pca Paramecium caudatum
# Ppr Paramecium primaurelia
# Pbi Paramecium biaurelia
# Pte Paramecium tetraurelia
# Ppe Paramecium pentaurelia
# Pse Paramecium sexaurelia

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'datapath=s' => \$dataPath,
			      'species=s' => \$species3abr,
			      'floating=i'  => \$floating));
die $usage if $help;
my $home = '/home/dsellis/';
$dataPath = $home.'data/IES_data/'; #default for local run
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
if($species3abr eq 'Pca' or
   $species3abr eq 'Pbi' or
   $species3abr eq 'Pte' or
   $species3abr eq 'Pse'){
    @lineage  = ('Eukaryota','Alveolata','Ciliophora','Intramacronucleata','Oligohymenophorea','Peniculida','Parameciidae','Paramecium');
}elsif($species3abr eq 'Tth'){
    @lineage  = ('Eukaryota','Alveolata','Ciliophora','Intramacronucleata', 'Oligohymenophorea', 'Hymenostomatida','Tetrahymenina', 'Tetrahymenidae', 'Tetrahymena');
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
    if($floating){
	$iesgffF = $dataPath.'internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.fl.gff3';
    }else{
	$iesgffF = $dataPath.'internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3';
    }
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
    if($floating){
	$iesgffF = $dataPath.'internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.fl.gff3';
    }else{
	$iesgffF = $dataPath.'internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3';
    }
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
    if($floating){
	$iesgffF = $dataPath.'internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3';
    }else{
	$iesgffF = $dataPath.'internal_eliminated_sequence_PGM_IES51.pt_51.gff3';
    }
}elsif($species3abr eq 'Pca'){
    $species = 'Paramecium caudatum';
    $taxonId = 5885;
    $speciesAbr = 'PCAU.43c3d.1.';
    $dataPath = $dataPath.'pcaudatum_43c3d_annotation_v2.0/';
    $cds = $dataPath.'pcaudatum_43c3d_annotation_v2.0.cds.fa';
    $protein = $dataPath.'pcaudatum_43c3d_annotation_v2.0.protein.fa';
    $gene = $dataPath.'pcaudatum_43c3d_annotation_v2.0.gene.fa';
    $gff3 = $dataPath.'pcaudatum_43c3d_annotation_v2.0.gff3';
    $scaffoldsF = $dataPath.'caudatum_43c3d_assembly_v1.fasta';
    $outputFile = $dataPath.'Pca.gnbk';
    if($floating){ # not yet available
	$iesgffF = $dataPath.'PCAUD_MIC10_IES.fl.gff3';
    }else{
	$iesgffF = $dataPath.'PCAUD_MIC10_IES.gff3';
    }
}elsif($species3abr eq 'Tth'){
    $species = 'Tetrahymena thermophila';
    $taxonId = 5911;
    $speciesAbr = 'TTHERM_';
    $dataPath = $dataPath.'tthermophila/';
    $cds = $dataPath.'T_thermophila_June2014_CDS.fasta';
    $protein = $dataPath.'T_thermophila_June2014_proteins.fasta';
    $gene = $dataPath.'T_thermophila_June2014_gene.fasta';
    $gff3 = $dataPath.'T_thermophila_June2014.gff3';
    $scaffoldsF = $dataPath.'T_thermophila_June2014_assembly.fasta';
    $outputFile = $dataPath.'Tth.gnbk';
    if($floating){ # not available
	$iesgffF = $dataPath.'';
    }else{
	$iesgffF = $dataPath.'';
    }
}else{
    die "unknown species";
}

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
	# if(defined($curGene)){
	# 	if(defined($CDSH{$curGene})){
	# 	    $entriesH{$prevScaffold}->add_SeqFeature($CDSH{$curGene});
	# 	    $curGene = $id[0]; #set the current gene
	# 	}else{
	# 	    $curGene = $id[0]; #for a gene without CDS
	# 	}
	# }else{
	# 	$curGene = $id[0];
	# }
	#build features
	#find all entries that are in the same scaffold
	$geneFeatureH{$id[0]} = new Bio::SeqFeature::Generic(-start       => $feature->start(),
							     -end         => $feature->end(),
							     -strand      => $feature->strand(),
							     -primary_tag => $feature->primary_tag(),
							     -tag => {gene     => $id[0]});
#	    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
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
#	    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$id[0]});
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

#add features to scaffolds
foreach my $scaffold(sort keys %scaffoldElementsH){
    foreach my $element (@{$scaffoldElementsH{$scaffold}}){
	if(defined($geneFeatureH{$element})){
	    $entriesH{$scaffold}->add_SeqFeature($geneFeatureH{$element});
#	    print $scaffold,' ';
#	    print $element,' ';
#	    print $geneFeatureH{$element}->start, ' ';
	    if(defined($CDSH{$element})){
#		print $CDSH{$element}->primary_tag();
		$entriesH{$scaffold}->add_SeqFeature($CDSH{$element});
	    }
#	    print "\n";
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
    if($species3abr eq 'Tth'){
	$proteinId = $geneIdNameH{$geneId};
    }elsif(($species3abr eq 'Pca') or
	   ($species3abr eq 'Pbi') or
	   ($species3abr eq 'Pte') or
	   ($species3abr eq 'Pse')){
	$proteinId = $geneId;
	$proteinId =~s/(.*)G(\d+)/$1P$2/ or die "$proteinId";
    }else{
	die;
    }
    return $proteinId;
}
