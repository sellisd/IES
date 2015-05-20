#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Getopt::Long;
my $help;
my $from;
my $to;
my $msa;
my $usage = <<HERE;

Change format of files using Bio::SeqIO
usage a2b.pl [OPTIONS] INPUTFILE OUTPUTFILE
where OPTIONS can be:
  -from:   input file format
  -to:     output file format
  -msa:    true if input is multiple sequence alignment (False by default)
  -help|?: this help screen
HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'from=s' => \$from,
			      'to=s'   => \$to,
			      'msa'    => \$msa,
		   ));
die $usage if $#ARGV < 1;
my $inputFile = $ARGV[0];
my $outputFile = $ARGV[1];

if(!$msa){
    my $in  = Bio::SeqIO->new(-file => $inputFile,
			      -format => $from);
    my $out = Bio::SeqIO->new(-file => '>'.$outputFile,
			      -format => $to);
    
    while ( my $seq = $in->next_seq() ) {
	$out->write_seq($seq);
    }
}else{
    my $in  = Bio::AlignIO->new(-file => $inputFile,
				-format => $from);
    my $out = Bio::AlignIO->new(-file => '>'.$outputFile,
				-format => $to);
    
    while ( my $aln = $in->next_aln() ) {
	$out->write_aln($aln);
    }
}

