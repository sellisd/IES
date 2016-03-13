#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
my $home = '/home/dsellis/';
#create multiple sequence alignments of silix groups
my $t_coffeeBin = '/panhome/sellis/tools/tcoffee/Version_11.00.8cbe486/bin/t_coffee';
#my $t_coffeeBin = 't_coffee';
my $dataPathNode = '/data/sellis/msas/';
my $dataPath = '/pandata/sellis/';
#my $dataPath = $home.'data/IES_data/';
my $fastaPath = $dataPath.'fasta/';
my $block = 10;

opendir(DH, $fastaPath) or die "$!: $fastaPath";
my @files = readdir(DH);
close DH;

#aligning sequences in parallel
print "running t_coffee for each group\n";

#run $block instances in each PBS file as each one is usualy too fast

my $counter = 0;
for(my $i = 0; $i <= $#files ; $i+=$block){
    open PBS, '>msa.'.$counter.'.pbs' or die $!;
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N msa.',$counter,"\n";
    print PBS '#PBS -e error.',$counter,"\n";
    print PBS '#PBS -o output.',$counter,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    print PBS 'date',"\n";
    for(my $j = 0; $j<$block; $j++){
	my $file = $files[$i+$j];
	next unless defined $file;
	next unless -f $fastaPath.$file;
	my $cmdl = $t_coffeeBin.' '.$fastaPath.$file.' -multi_core no  &> '.$file.'.log';
 	print PBS "$cmdl","\n";    
    }
    print PBS "mkdir $dataPathNode\n"; #make sure the file exists
    print PBS "mv *.log *.aln *.dnd *.html $dataPathNode\n";
    print PBS 'echo telos',"\n";
    print PBS 'date',"\n";
    close PBS;
    system "qsub msa.".$counter.".pbs";
#    print "qsub msa.".$counter.".pbs\n";
    $counter++;
 #   my $waiting = `qstat |grep q1hour|awk '\$5!="C"'|wc`;
 #   while((split " ", $waiting)[0] > 2900){ #buffer of 100
 #	$waiting = `qstat |grep q1hour|awk '\$5!="C"'|wc`;
 #	print "waiting for queue to clean up $waiting\n";
 #	sleep 1;
 #   }
}


