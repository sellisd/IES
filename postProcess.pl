#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
use Getopt::Long;
use Bio::Species;
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

if($speciesAbr eq 'Ppr'){
    $speciesName = 'Paramecium primaurelia';
    $taxonId = 5886;
}elsif($speciesAbr eq 'Pbi'){
    $speciesName = 'Paramecium biaurelia';
    $taxonId = 65126;
#Paremecium biaurelia
}elsif($speciesAbr eq 'Pte'){
    $speciesName = 'Paramecium tetraurelia';
    $taxonId = 5888;
#Paramecium tetraurelia
}elsif($speciesAbr eq 'Pen'){
    $speciesName = 'Paramecium pentaurelia';
    $taxonId = 43138;
#Paramecium pentaurelia
}elsif($speciesAbr eq 'Pse'){
    $speciesName = 'Paramecium sexaurelia';
    $taxonId = 65128;
#Paramecium sexaurelia
}elsif($speciesAbr eq 'Pca'){
    $speciesName = 'Paramecium caudatum';
    $taxonId = 5885;
}elsif($speciesAbr eq 'Tth'){
    $speciesName = 'Tetrahymena thermophila';
    $taxonId = 5911;
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
