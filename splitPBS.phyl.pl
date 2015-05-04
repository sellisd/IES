#!/usr/bin/perl
use warnings;
use strict;
#if a pbs is too long for execution in the cluster split the file to smaller pieces and rerun it

# find which ones have non zero error file
# read their corresponding pbs file
# find the listfile and split it in pieces, for each piece make a new pbs file

my $path = $ARGV[0];
opendir(DH, $path) or die $!;
my @files = readdir(DH);
close DH;
my $step = 1;
my @files2sub;
foreach my $file (@files){
    next unless -s $path.$file; #only error files with not nzero size
    next unless $file =~ /error.(\d+)/; #keep only error files
    my $number = $1;
    open PBS, $path.'phyl.'.$number.'.pbs' or die $!;
    my @head;
    my @commands;
    my @tail;
    while(my $line = <PBS>){
	if (substr($line,0,1) eq '#'){
#if q1hour make q1week
	    chomp $line;
#	    if ($line eq '#PBS -q q1hour'){
#		$line = '#PBS -q q1week';
#	    }
	    push @head,$line."\n";
	}elsif(substr($line,0,1) eq '/'){
	    die unless $line =~ /listfiles\.(\d+)/; # find listfile and split it in smaller chuncks
	    my $number = $1;
	    open LF, '/pandata/sellis/msas/listfiles/listfiles.'.$number or die $!;
	    my $lineCounter = 0;
	    my $newlistfile = '/pandata/sellis/msas/listfiles/listfiles.'.$number.'.'.$lineCounter;
	    while (my $l = <LF>){
		open LFN, '>'.$newlistfile or die $!;
		print LFN $l;
		close LFN;
		#make one PBS file for each chunk	   
		$line = '/pandata/penel/exemple_raxml/ms_fasta_files_2_srs.pl  --option=\'-m GTRGAMMA\' --type=dna -modele=DNA -partition /pandata/sellis/msas/listfiles/listfiles.'.$number.'.'.$lineCounter.' /pandata/sellis/msas/trees/out.'.$number.'.'.$lineCounter.'.srs'."\n";
		#foreach piece make a new pbs file:
		my $fileName = "phyl.".$number.'.'.$lineCounter.".pbs";
		open NBPS, '>'.$fileName or die $!;
		foreach my $header (@head){
		    if($header =~/#PBS(.+)\.$number/){
			chomp $header;
			print NBPS $header.'.'.$lineCounter."\n";
		    }else{
			print NBPS $header;
		    }
		}
		print NBPS $line;
		close NBPS;
		print $fileName,"\n";
		$lineCounter++;
#		system "qsub $sub";
	    }
	    close LF;
	}
    }
}
