#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

#create multiple sequence alignments of silix groups
my $t_coffeeBin = '/panhome/sellis/tools/tcoffee/Version_11.00.8cbe486/bin/t_coffee';
my $dataPath = '/pandata/sellis/';
my $fastaPath = '/pandata/sellis/msas/fasta/';
my $logfilePath = '/pandata/sellis/msas/log/';
my $alnfilePath = '/pandata/sellis/msas/aln/';
my $dndfilePath = '/pandata/sellis/msas/dnd/';
my $htmlfilePath = '/pandata/sellis/msas/html/';

#make required directories
mkdir $logfilePath unless (-d $logfilePath);
mkdir $dndfilePath unless (-d $dndfilePath);
mkdir $htmlfilePath unless (-d $htmlfilePath);
mkdir $alnfilePath unless (-d $alnfilePath);

opendir(DH, $fastaPath) or die $!;
my @files = readdir(DH);
close DH;

#aligning sequences in parallel
print "running t_coffee for each group\n";

my $counter = 0;
foreach my $file (@files){
    next unless -f $fastaPath.$file;
    #use single core for each group
    my $cmdl = $t_coffeeBin.' '.$fastaPath.$file.' -multi_core no  &> '.$file.'.log';
#make PBS file
    open PBS, '>'.$file.'.pbs' or die $!;
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N msa.',$file,"\n";
    print PBS '#PBS -e error.',$file,"\n";
    print PBS '#PBS -o output.',$file,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    print PBS "$cmdl","\n";
    print PBS "mv *.log $logfilePath\n";
    print PBS "mv *.aln $alnfilePath\n";
    print PBS "mv *.dnd $dndfilePath\n";
    print PBS "mv *.html $htmlfilePath\n";
    print PBS 'echo telos',"\n";
    close PBS;
#    $pm->finish;
    system "qsub ".$file.".pbs";
    my $waiting = 3000; #initialize large
    while((split " ", $waiting)[0] < 2900){ #buffer of 100
	$waiting = `qstat |grep q1hour|awk '$5!="C"'|wc`;
	print "waiting for queue to clean up $waiting\n";
	sleep 1;
    }
    $counter++; #make sure there is enough space to submit new jobs
}
#$pm->wait_all_children;

