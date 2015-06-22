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

Change format of files using Bio::SeqIO or Bio::AlignIO from input file format to output.
usage a2b.pl [OPTIONS] INPUTFILE(S)
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
die $usage if $#ARGV < 0;
die $usage if $help;
foreach my $inputFile (@ARGV){
    my $outputFile = $inputFile.'.'.$to; 
    die "$outputFile" if -f $outputFile; #do not overwrite any files
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
}
