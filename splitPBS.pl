#!/usr/bin/perl
use warnings;
use strict;
#if a pbs is too long for execution in the cluster split the file to smaller pieces and rerun it
#find which ones to split
opendir(DH, $ARGV[0]) or die $!;
my @files = readdir(DH);
close DH;
my $step = 10;
foreach my $file(@files){
    next unless -s $file;
    next unless $file =~ /error.(\d+)/; #keep only error files
    print $file;
    my $number = $1;
    open PBS, 'msa.'.$number.'.pbs' or die $!;
    my @head;
    my @commands;
    my @tail;
    while(my $line = <PBS>){
	if (substr($line,0,1) eq '#'){
	    push @head,$line;
	}elsif(substr($line,0,1) eq '/'){
	    push @commands,$line;
	}else{
	    push @tail, $line;
	}
    }
    close PBS;
#split file in pieces
    my $counter = 0;
    for(my $i = 0; $i<$#commands; $i+=$step){
	my $fileName = "msa.".$number.'.'.$counter.".pbs";
	print $fileName,"\n";
	open NBPS, '>'.$fileName or die $!;
	foreach my $header (@head){
	    if($header =~/#PBS(.+)\.$number/){
		$header = $header.'.'.$counter;
	    }
	    print NBPS $header;
	}
	my $to = ($i+$step>$#commands?$#commands:$i+$step);
	my @slice = @commands[$i..$to];
	print NBPS @slice;
	print NBPS @tail;
	$counter++;
	close NBPS;
	print $fileName,"\n";
#	system "qsub $fileName";
    }
}
