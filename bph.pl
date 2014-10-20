#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Parallel::ForkManager;

#script runs all-against-all blast of proteins after filtering low complexity regions and performs silix grouping

#make sure each protein has unique ID including in which species it belongs
#concatenate all fasta files to create on big one
#cat pbiaurelia/biaurelia_V1-4_annotation_v1.protein.fa psexaurelia/sexaurelia_AZ8-4_annotation_v1.protein.fa ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.4.protein.fa > combined.fa

#make database for all proteins
#  makeblastdb -in pb.fa -dbtype prot -out combined

# for some reason cat messed up some new lines so watch the output of makeblsatdb

my $cores  = 8;
my $pm = Parallel::ForkManager->new($cores);
# pipeline for analysis of Paremecium data

#set data paths
my $blastdbPath = 'working/';
my $fastaInPath = 'working/';
my $fastaOutPath = 'working/';
my $blastOutPath = 'working/';

my @tempFiles; #fill array with temporary files to be deleted

#blast all against all but keep all hits and output in tab file
#  blastp -query combined.fa -db combined.fa -outfmt 6
#or ?faster? read fasta file one sequence at a time and blast it against the database

my $file = "/Users/diamantis/data/IES_data/working/combined.fa";
my $seqF = Bio::SeqIO->new('-file'=>$file,
			   '-format'=>'fasta');
#split fasta file in chunkcs
my $batch = 100;
my $counter = 0;
my $fileCounter = 0;
my $fastaOutFile = $fastaOutPath.'chunk.'.$fileCounter.'.fa';
my $seqOut = Bio::SeqIO->new('-file'=>'>'.$fastaOutFile,
			     '-format'=>'fasta');
push @tempFiles, $seqOut;
print "splitting fasta file in chunks of $batch sequences\n";
while(my $seqO = $seqF->next_seq){
    print "   $counter\r";
    if ($counter % $batch == 0){
	$fastaOutFile = $fastaOutPath.'chunk.'.$fileCounter.'.fa';
	$seqOut = Bio::SeqIO->new('-file'=>'>'.$fastaOutFile,
				  '-format'=>'fasta');
	push @tempFiles, $fastaOutFile;
	$fileCounter++;
    }
    $seqOut->write_seq($seqO);
    $counter++;
}
print "   $counter sequences in $fileCounter pieces\n";
print "     done\n";

print "running blast\n";
for(my $i = 0 ; $i<$fileCounter; $i++){
    print "   $i/$fileCounter\r";
    my $pid = $pm->start and next;
    my $cmdString = 'blastp -query '.$fastaOutPath.'chunk.'.$i.'.fa -db '.$fastaInPath.'combined -outfmt 6';
#    my $cmdString = '/home/dsellis/tools/ncbi-blast-2.2.29+/bin/blastp -query '.$fastaOutPath.'chunk.'.$i.'.fa -db '.$fastaInPath.'combined -outfmt 6';
    my $result = `$cmdString`;
    my $blastOutFile = $blastOutPath.'blastout.'.$i.'.dat';
    open OUT, '>'.$blastOutFile or die $!;
    print OUT $result;
    close OUT;
    push @tempFiles, $blastOutFile;
    $pm->finish;
}
$pm->wait_all_children;
print "     done\n";

print "stitch output\n";
my $cmdString = 'cat '.$blastOutPath.'blastout.*.dat > '.$blastOutPath.'4silix.dat';
system($cmdString);
#blast one chunk at a time (on a different core 8) ) 
#output in tabular format for Silix

# run silix for clustering
#silix fastafile blastfile > silix.output
#silix working/combined.fa working/4silix.dat > working/silix.output
#
#delete intermidiate files
print "deleting temp files\n";
unlink @tempFiles or die $!;

