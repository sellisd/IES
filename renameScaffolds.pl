#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
#read a genbank file and rename the scaffold (accession and locus) to be compatible with the acnuc database
my $inputfile =$ARGV[0];
my $outputfile=$ARGV[1];
my $species3abr = $ARGV[2];

my $dataIn = Bio::SeqIO->new('-file' => $inputfile,
			      '-format' => 'genbank');
my $dataOut = Bio::SeqIO->new('-file' => '>'.$outputfile,
			       '-format' => 'genbank');

#count scaffolds
my $NumberOfScaffolds = 0;
print "counting scaffolds\n";
while(my $scaffoldSeq = $dataIn->next_seq){
    $NumberOfScaffolds++;
}
#re-read (not elegant)
$dataIn = Bio::SeqIO->new('-file' => $inputfile,
			  '-format' => 'genbank');
print "renaming\n";
#rename scaffolds with acnuc compatible names
my $NumberOfDigits = length($NumberOfScaffolds + 1);
while(my $scaffoldSeq = $dataIn->next_seq){
    my $scname = $scaffoldSeq->accession_number();
    $scname =~ /sca?ff?o?l?d?.*?_(\d+)/;
    my $number = $1;
    my $printString = '%0'.$NumberOfDigits.'d';
    my $Padded = $species3abr.'_'.sprintf($printString,$number);
    $scaffoldSeq->display_id($Padded);
    $scaffoldSeq->accession_number($Padded);
    $dataOut->write_seq($scaffoldSeq);
}

