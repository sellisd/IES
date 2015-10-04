#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Getopt::Long;

my $help;
my $input;
my $output;
my $addTA;
my $usage = <<HERE;

Create a fasta file with IES sequence from gff3 files.
usage:
  
ies2fasta.pl [OPTIONS] IN OUT
where OPTIONS can be:
    -addTA:  add a trailing TA to all sequences (for P. caudatum IES)
    -help|?: this help screen
and
   IN: input file
  OUT: output file
 
HERE

die $usage unless (GetOptions('help|?' => \$help,
		   'addTA' => \$addTA));
die $usage if $#ARGV < 1;
die $usage if $help;

#read ies gff
$input = $ARGV[0]; #'/home/dsellis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51_with_ies.gff3';
$output = $ARGV[1]; #'>/home/dsellis/data/IES_data/ptetraurelia/ies.fa';

if(! -e $input){
    die "file does not exist: $input\n";
}
if(-e $output){
    die "output file already exists, didn't overwrite: $output\n";
}

my $outS = Bio::SeqIO->new('-file' => '>'.$output,
			   '-format' => 'fasta');
my $inS = Bio::Tools::GFF->new('-file' => $input,
			       '-gff_version'=> 3);
while(my $feature = $inS->next_feature()){
    my $id = ($feature->get_tag_values('ID'))[0];
    my $seq =  ($feature->get_tag_values('sequence'))[0];
    if ($addTA){
	$seq .= 'TA';
    }
    my $seqO = Bio::Seq->new('-seq'      => $seq,
			     '-id'       => $id,
			     '-alphabet' => 'dna');
    $outS->write_seq($seqO);			 
}
