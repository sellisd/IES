#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
use Getopt::Long;
my $help;
my $usage = <<HERE;

Add IES information in a genbank file
usage: postProcess.pl IES.gff3 input.gnbk

HERE
die $usage unless (GetOptions('help|?' => \$help));
die $usage if $help;
my $iesgffF = $ARGV[0];#'/Users/diamantis/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3';
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
      my $newFeature = new Bio::SeqFeature::Generic(-start => $ies->{'start'},
						    -end   => $ies->{'end'},
						    -primary_tag   => "IES_junction",
						    -tag   => {'score' => $ies->{'score'},
							       'sequence' => $ies->{'sequence'},
							       'id'  => $ies->{'id'}
						    }
	  );
      $seqO->add_SeqFeature($newFeature);
      my $iesSeqO = Bio::Seq->new('-display_id' => $ies->{'id'},
				  '-format' =>'genbank',
				  '-accession_number' => $ies->{'id'},
				  '-start' => $ies->{'start'},
				  '-end' => $ies->{'end'},
				  '-alphabet' => 'dna',
				  '-accession_number' => $ies->{'id'},
		#		  '-species' => 'not implememted yet!',
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
#TODO add species and make ACCESSION valid (check also restrictions for LOCUS)
#make a string with species abr.scaffold and incrementing number
      $iesgnbkOut->write_seq($iesSeqO);
    }
    $gnbkOut->write_seq($seqO);

}
