#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
#use Parallel::ForkManager; not useful in the cluster

#script runs all-against-all blast of proteins after filtering low complexity regions and performs silix grouping


my $file = "/pandata/sellis/combined.fa"; #fasta file
#my $file = "/pandata/sellis/test.fa"; #fasta file
my $blastpBin = '/usr/remote/ncbi-blast-2.2.29+/bin/blastp';
my $silixBin = 'silix/usr/remote/silix-1.2.8/bin/';

#my $cores  = 48;
#my $pm = Parallel::ForkManager->new($cores);
# pipeline for analysis of Paremecium data

#set data paths
my $blastdbPath = '/pandata/sellis/';
my $fastaOutPath = '/pandata/sellis/';
my $blastOutPath = '/pandata/sellis/';
my $silixOutPath = '/pandata/sellis/';
my @tempFiles; #fill array with temporary files to be deleted

#blast all against all but keep all hits and output in tab file
#  blastp -query combined.fa -db combined.fa -outfmt 6
#or ?faster? read fasta file one sequence at a time and blast it against the database


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
#    my $pid = $pm->start and next;
#write pbs file
    open PBS, '>'.$i.'.pbs' or die $!;
    my $cmdString = $blastpBin.' -query '.$fastaOutPath.'chunk.'.$i.'.fa -db '.$blastdbPath.'combined -outfmt 6 -out '.$blastOutPath.'blastout.'.$i.'.dat -seg yes';
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N blastp.',$i,"\n";
    print PBS '#PBS -e error.',$i,"\n";
    print PBS '#PBS -o output.',$i,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    print PBS "$cmdString","\n";
    print PBS 'echo telos',"\n";
    close PBS;

#submit file
   # my $result = `$cmdString`;
    print "qsub ".$i.'.pbs',"\n";
    system "qsub ".$i.'.pbs';
   # my $blastOutFile = $blastOutPath.'blastout.'.$i.'.dat';
   # open OUT, '>'.$blastOutFile or die $!;
   # print OUT $result;
   # close OUT;
#    $pm->finish;
}
#$pm->wait_all_children;
print "     done\n";



