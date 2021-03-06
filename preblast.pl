#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use File::Spec::Functions qw(catfile);
my $file = $ARGV[0];
my $fastaOutPath = $ARGV[1];

mkdir $fastaOutPath unless -d $fastaOutPath;
#prepare fasta files for running blast on cluster
my $seqF = Bio::SeqIO->new('-file'=>$file,
			   '-format'=>'fasta');
#split fasta file in chunkcs
my $batch = 200;
my $counter = 0;
my $fileCounter = 0;
my $fastaOutFile = $fastaOutPath.'chunk.'.$fileCounter.'.fa';
my $seqOut = Bio::SeqIO->new('-file'=>'>'.$fastaOutFile,
			     '-format'=>'fasta');

print "splitting fasta file in chunks of $batch sequences\n";
while(my $seqO = $seqF->next_seq){
    print "   $counter\r";
    if ($counter % $batch == 0){
	$fastaOutFile = catfile($fastaOutPath,'chunk.'.$fileCounter.'.fa');
	$seqOut = Bio::SeqIO->new('-file'=>'>'.$fastaOutFile,
				  '-format'=>'fasta');
	$fileCounter++;
    }
    $seqOut->write_seq($seqO);
    $counter++;
}

print "   $counter sequences in $fileCounter pieces\n";
print "     done\n";
