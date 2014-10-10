#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
# pipeline for analysis of Paremecium data

#make sure each protein has unique ID including in which species it belongs
#concatenate all fasta files to create on big one
#cat p1.fa s2.fa > combined.fa

#make database for all proteins
#  makeblastdb -in combined.fa -dbtype prot

# for some reason cat messed up some new lines so watch the output of makeblsatdb

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
my $seqOut = Bio::SeqIO->new('-file'=>'>chunk.'.$fileCounter,
			     '-format'=>'fasta');
while(my $seqO = $seqF->next_seq){
  if ($counter % $batch == 0){
    $seqOut = Bio::SeqIO->new('-file'=>'>chunk.'.$fileCounter,
			     '-format'=>'fasta');
    $fileCounter++;
  }
  $seqOut->write_seq($seqO);
  $counter++;
}

#blast one chunk at a time (on a different core 8) ) 

#output in tabular format for Silix

# run silix for clustering
#silix fastafile blastfile > silix.output

