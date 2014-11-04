#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

#create multiple sequence alignments of silix groups
my $t_coffeeBin = '/panhome/sellis/tools/tcoffee/Version_11.00.8cbe486/bin/t_coffee';
#my $t_coffeeBin = 't_coffee';
my $dataPath = '/pandata/sellis/';
#my $dataPath = '/Users/diamantis/data/IES_data/';
my $fastaPath = $dataPath.'msas/fasta/';
my $logfilePath = $dataPath.'msas/log/';
my $alnfilePath = $dataPath.'msas/aln/';
my $dndfilePath = $dataPath.'msas/dnd/';
my $htmlfilePath = $dataPath.'msas/html/';
my $block = 500;
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

#run $block instances in each PBS file as each one is usualy too fast
my $counter = 0;
for(my $i = 0; $i < $#files ; $i+=$block){
    open PBS, '>msa.'.$counter.'.pbs' or die $!;
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N msa.',$counter,"\n";
    print PBS '#PBS -e error.',$counter,"\n";
    print PBS '#PBS -o output.',$counter,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    for(my $j = 0; $j<$block; $j++){
	my $file = $files[$i+$j];
	next unless defined $file;
	next unless -f $fastaPath.$file;
	my $cmdl = $t_coffeeBin.' '.$fastaPath.$file.' -multi_core no  &> '.$file.'.log';
 	print PBS "$cmdl","\n";    
    }
    print PBS "mv *.log $logfilePath\n";
    print PBS "mv *.aln $alnfilePath\n";
    print PBS "mv *.dnd $dndfilePath\n";
    print PBS "mv *.html $htmlfilePath\n";
    print PBS 'echo telos',"\n";
    close PBS;
    $counter++;
    system "qsub msa.".$counter.".pbs";
    my $waiting = `qstat |grep q1hour|awk '\$5!="C"'|wc`;
    while((split " ", $waiting)[0] > 2900){ #buffer of 100
 	$waiting = `qstat |grep q1hour|awk '\$5!="C"'|wc`;
 	print "waiting for queue to clean up $waiting\n";
 	sleep 1;
    }
}


