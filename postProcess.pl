#!/usr/bin/perl
use warnings;
use strict;
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
die "@annotations" unless defined($id);
  $IESH{$id} = {'scaffold' => $scaffold,
		'start' => $start,
		'end'   => $end,
		'score' => $score,
		'sequence' => $sequence};
}
close IN;
print Dumper %IESH;
