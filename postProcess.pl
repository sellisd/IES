#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
# add IES in the constructed genbank file
my $iesgffF = '/Users/diamantis/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3';
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
	       'sequence' => $sequence};
  if(defined($IESH{$scaffold})){
      push @{$IESH{$scaffold}}, $entry;
  }else{
      $IESH{$scaffold} = [$entry];
  }
}
close IN;

#read genbank file
#go through sequencies and add as features the IES

my $genbank = 'output.gnbk';
my $gnbkIn = Bio::SeqIO->new('-file' => $genbank,
 			    '-format' => 'genbank');

my $gnbkOut = Bio::SeqIO->new('-file'=> '>newOut.gnbk',
			      '-format' => 'genbank');
while(my $seqO = $gnbkIn -> next_seq()){
  my $scaffold = $seqO->display_id();
#loop through iess in this scaffold
  foreach my $ies (@{$IESH{$scaffold}}){
      my $newFeature = new Bio::SeqFeature::Generic(-start => $ies->{'start'},
						    -end   => $ies->{'end'},
						    -primary_tag   => "IES_junction",
						    -tag   => {'score' => $ies->{'score'},
							       'sequence' => $ies->{'sequence'}
						    }
	  );
      $seqO->add_SeqFeature($newFeature);
  }
  $gnbkOut->write_seq($seqO);
  #find which ones are in this scaffold
#make features
#add feature to scaffold
}
