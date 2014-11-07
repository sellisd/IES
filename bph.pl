#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

#script runs all-against-all blast of proteins after filtering low complexity regions and performs silix grouping
my $blastpBin = '/usr/remote/ncbi-blast-2.2.29+/bin/blastp';
my $dataPath = '/pandata/sellis/working/fasta/'; #@pbil
my $nodeDataPath = '/data/sellis/working/fasta/'; #@nodes
my $blastdbPath = '/data/sellis/working/';
my $blastOutPath = '/data/sellis/working/blastout/';

mkdir $blastOutPath unless -d $blastOutPath;

opendir(DH, $dataPath) or die $!;
my @files = readdir(DH);
close DH;
my $fileCounter = 0;
foreach my $file (@files){
    next unless -f $dataPath.$file;
    print $file,"\n";
    open PBS, '>'.$fileCounter.'.pbs' or die $!;
    my $cmdString = $blastpBin.' -query '.$nodeDataPath.$file.' -db '.$blastdbPath.'combined -outfmt 6 -out '.$blastOutPath.'blastout.'.$file.'.dat -seg yes';
    #write pbs file
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N blastp.',$fileCounter,"\n";
    print PBS '#PBS -e error.',$fileCounter,"\n";
    print PBS '#PBS -o output.',$fileCounter,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    print PBS "$cmdString","\n";
    print PBS 'echo telos',"\n";
    close PBS;
    print "qsub ".$fileCounter.'.pbs',"\n";
    system "qsub ".$fileCounter.'.pbs';
    $fileCounter++;
}
