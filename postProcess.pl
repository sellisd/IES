#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
use Getopt::Long;
use Bio::Species;
use lib'.';
use functions;
my $help;
my $speciesAbr;
my $usage = <<HERE;

Add IES information in a genbank file
  usage: postProcess.pl -species Pab IES.gff3 input.gnbk

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'species=s' => \$speciesAbr));
die $usage if $help;
my $taxonId;
my $speciesName;
my @lineage = ('Eukaryota','Alveolata','Ciliophora','Intramacronucleata','Oligohymenophorea','Peniculida','Parameciidae','Paramecium');

my $addTA = 1; #always true in new dataset

if($speciesAbr eq 'ppr'){
    $speciesName = 'Paramecium primaurelia';
}elsif($speciesAbr eq 'pbi'){
    $speciesName = 'Paramecium biaurelia';
}elsif($speciesAbr eq 'pte'){
    $speciesName = 'Paramecium tetraurelia';
}elsif($speciesAbr eq 'ppe'){
    $speciesName = 'Paramecium pentaurelia';
}elsif($speciesAbr eq 'pse'){
    $speciesName = 'Paramecium sexaurelia';
}elsif($speciesAbr eq 'poc'){
    $speciesName = 'Paramecium octaurelia';
}elsif($speciesAbr eq 'ptr'){
    $speciesName = 'Paramecium tredecaurelia';
}elsif($speciesAbr eq 'pso'){
    $speciesName = 'Paramecium sonneborni';
}elsif($speciesAbr eq 'pca'){
    $speciesName = 'Paramecium caudatum';
}elsif($speciesAbr eq 'tth'){
    $speciesName = 'Tetrahymena thermophila';
}else{
    print 'Not known species abreviation: $speciesAbr',"\n";
    die $usage;
}

push @lineage, $speciesName;
@lineage = reverse(@lineage);
my $speciesO = Bio::Species->new(-classification => \@lineage);

my $iesgffF = $ARGV[0];
my $genbank = $ARGV[1]; 
my $genBankOut = $genbank;
$genBankOut =~ s/\.gnbk/.IES.gnbk/;
my $iesgnbkOutF = $genBankOut;
$iesgnbkOutF =~ s/.IES.gnbk/.ies/;
my $IESgff = Bio::Tools::GFF->new('-file' => $iesgffF,
				  '-format' => 'gff3');
# parse gff3 with annotations
my %IESH;
open IN, $iesgffF;
while (my $line = <IN>){
    next if substr($line, 0, 1) eq '#'; #skip comments
    chomp $line;
    my @ar = split "\t", $line;
    my $scaffold = $ar[0];
    my $start = $ar[3];
    my $end = $ar[4];
    my $score = $ar[5];
    my $annotation = $ar[8];
    my $id;
    my $sequence;
    my $alt_seq;
    my @annotations = split ';', $annotation;
    my $counter = 0;
    foreach my $annot (@annotations){
	if ($annot =~ /ID=(.*)/){
	    $id = $1;
	}elsif($annot =~ /sequence=(.*)/){
	    $sequence = $1;
	    if($addTA){
		$sequence = $sequence.'TA';
	    }
	}elsif($annot =~ /alternative_IES_seq=(.*)/){
	    $alt_seq = $1;
	}
	$counter++;
    }
    my $entry = {'id' => $id,
		 'start' => $start,
		 'end'   => $end,
		 'score' => $score,
		 'sequence' => $sequence,
		 'alt_sequence' => $alt_seq};
    if(defined($IESH{$scaffold})){
	push @{$IESH{$scaffold}}, $entry;
    }else{
	$IESH{$scaffold} = [$entry];
    }
}
close IN;

#read genbank file
#go through sequencies and add as features the IES

my $gnbkIn = Bio::SeqIO->new('-file' => $genbank,
    '-format' => 'genbank');

my $gnbkOut = Bio::SeqIO->new('-file'=> '>'.$genBankOut,
    '-format' => 'genbank');
my $iesgnbkOut = Bio::SeqIO->new('-file' => '>'.$iesgnbkOutF,
    '-format'=> 'genbank');

while(my $seqO = $gnbkIn -> next_seq()){
    my $scaffold = $seqO->display_id();
#loop through iess in this scaffold
    foreach my $ies (@{$IESH{$scaffold}}){
	$scaffold =~ /.*?_(\d+)/;
	my $scaffoldNo = $1 or die $scaffold;
	my $accNumber = $speciesAbr.$scaffoldNo.'_'.$ies->{'start'};
	
	my $newFeature = new Bio::SeqFeature::Generic(-start => $ies->{'start'},
						      -end   => $ies->{'end'},
						      -primary_tag   => "IES_junction",
						      -tag   => {'score' => $ies->{'score'},
#								 'sequence' => $ies->{'sequence'}, #do not include sequence, it will be part of a separate file
								 'id'  => $ies->{'id'},
								 'db_xref' => $accNumber
						      }
	    );
	$seqO->add_SeqFeature($newFeature);
	my $iesSeqO = Bio::Seq->new('-display_id' => $accNumber,
				    '-desc' => "$speciesName IES located on scaffold $scaffoldNo at $ies->{'start'} (MIC coordinates)",
				    '-format' =>'genbank',
				    '-start' => $ies->{'start'},
				    '-end' => $ies->{'end'},
				    '-alphabet' => 'dna',
				    '-accession_number' => $accNumber,
				    '-species' => $speciesO,
				    '-seq' => $ies->{'sequence'}
	    );
	my $iesFeature = new Bio::SeqFeature::Generic(-start => $ies->{'start'},
						      -end   => $ies->{'end'},
						      -primary_tag => 'IES',
						      -tag   => {'score' => $ies->{'score'},
								 'id'  => $ies->{'id'},
								 'scaffold' => $scaffold}
	    );
	if(defined($ies->{'alt_sequence'})){
	    $iesFeature->add_tag_value('alternative sequence',$ies->{'alt_sequence'});
	}
	$iesSeqO->add_SeqFeature($iesFeature);
#make a string with species abr.scaffold and incrementing number
#add definition?
	$iesgnbkOut->write_seq($iesSeqO);
    }
    $gnbkOut->write_seq($seqO);
    
}
