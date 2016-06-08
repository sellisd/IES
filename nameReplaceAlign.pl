#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use lib'.';
use functions;
use Getopt::Long;
my $help;
my $usage = <<HERE;

Replace gene names with species names in multiple alignment files
usage nameReplaceAlign.pl [OPTIONS] INPUTFILE(S)
where OPTIONS can be:
  -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?'  => \$help));
die $usage if $#ARGV < 0;
die $usage if $help;

my @species = qw/
Paramecium_primaurelia
Paramecium_biaurelia
Paramecium_tetraurelia
Paramecium_pentaurelia
Paramecium_sexaurelia
Paramecium_octaurelia
Paramecium_tredecaurelia
Paramecium_sonneborni
Paramecium_caudatum
Tetrahymena_thermophila
/;

foreach my $inputFile (@ARGV){
    print $inputFile,"\n";
    my $in = Bio::SeqIO->new(-file   => $inputFile,
			     -format => 'Fasta');
    my $outputFile = $inputFile.'.renamed';
    my $out = Bio::SeqIO->new(-file => '>'.$outputFile,
			      -format => 'Fasta');
    my $hasTth = 0;
    while(my $seqO = $in->next_seq()){
	my $species = gene2species($seqO->id);
	$hasTth = 1 if $species eq 'Tetrahymena_thermophila';
	my $seq = $seqO->seq();
	my $newSeqO = Bio::Seq->new(-id => $species,
				    -seq => $seq);
	$out->write_seq($newSeqO);
    }
}
