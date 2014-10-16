#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Parallel::ForkManager;
#make sure each protein has unique ID including in which species it belongs
#concatenate all fasta files to create on big one
#cat p1.fa s2.fa > combined.fa

#make database for all proteins
#  makeblastdb -in combined.fa -dbtype prot

# for some reason cat messed up some new lines so watch the output of makeblsatdb

my $pm = Parallel::ForkManager->new(8);
# pipeline for analysis of Paremecium data

#set data paths
my $blastdbPath = 'working/';
my $fastaInPath = 'working/';
my $fastaOutPath = 'working/';
my $blastOutPath = 'working/';

#blast all against all but keep all hits and output in tab file
#  blastp -query combined.fa -db combined.fa -outfmt 6
#or ?faster? read fasta file one sequence at a time and blast it against the database
my $file = "/Users/diamantis/data/IES_data/working/combined.fa";
my $seqF = Bio::SeqIO->new('-file'=>$file,
			   '-format'=>'fasta');
#split fasta file in chunkcs
my $batch = 10;
my $counter = 0;
my $fileCounter = 0;
my $seqOut = Bio::SeqIO->new('-file'=>'>'.$fastaOutPath.'chunk.'.$fileCounter.'.fa',
			     '-format'=>'fasta');
print "splitting fasta file in chunks of $batch sequences\n";
while(my $seqO = $seqF->next_seq){
    print "   $counter\r";
    if ($counter % $batch == 0){
	$seqOut = Bio::SeqIO->new('-file'=>'>'.$fastaOutPath.'chunk.'.$fileCounter.'.fa',
				  '-format'=>'fasta');
	$fileCounter++;
    }
    $seqOut->write_seq($seqO);
    $counter++;
}
print "   $counter sequences in $fileCounter pieces\n";
print "   done\n";

print "running blast\n";
for(my $i = 0 ; $i<$fileCounter; $i++){
    print "   $i/$fileCounter\r";
    my $pid = $pm->start and next;
    my $cmdString = 'blastp -query '.$fastaOutPath.'chunk.'.$i.'.fa -db '.$fastaInPath.'combined.fa -outfmt 6';
    my $result = `$cmdString`;
    open OUT, '>'.$blastOutPath.'blastout.'.$i.'.dat' or die $!;
    print OUT $result;
    close OUT;
    $pm->finish;
}
print "   done\n";

print "stitch output\n";
my $cmdString = 'cat '.$blastOutPath.'blastout.*.dat > '.$blastOutPath.'4silix.dat';
system($cmdString);
#blast one chunk at a time (on a different core 8) ) 
#output in tabular format for Silix

# run silix for clustering
#silix fastafile blastfile > silix.output

#delete intermidiate files
