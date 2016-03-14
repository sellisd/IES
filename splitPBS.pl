#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);

#if a pbs is too long for execution in the cluster split the file to smaller pieces and rerun it
#find which ones to split
my $path = $ARGV[0];
opendir(DH, $path) or die $!;
my @files = grep{/error\.\d+/} readdir(DH);
close DH;
my $step = 1;
my @files2sub;
foreach my $file (@files){
    next unless -s catfile($path, $file); #only error files with not zero size
    $file =~ /error.(\d+)/; # extract number
    my $number = $1;
    #if error file an expected warning
    open ERR, catfile($path, $file) or die $!;
    my $line = readline(ERR);
    chomp $line;
    next if $line eq 'mkdir: cannot create directory ‘/data/sellis/msas/’: File exists';
    close ERR;
    open PBS, catfile($path,'msa.'.$number.'.pbs') or die $!;
    my @head;
    my @commands;
    my @tail;
    while(my $line = <PBS>){
	if (substr($line,0,1) eq '#'){
#if q1hour make q1week
	    chomp $line;
	    if ($line eq '#PBS -q q1hour'){
#		$line = '#PBS -q q1week';
	    }
	    push @head,$line."\n";
	}elsif(substr($line,0,1) eq '/'){
	    push @commands,$line;
	}else{
	    push @tail, $line;
	}
    }
    close PBS;
#split file in pieces
    my $counter = 0;
    for(my $i = 0; $i<=$#commands; $i+=$step){
	my $fileName = catfile($path, 'msa.'.$number.'.'.$counter.'.pbs');
	open NBPS, '>'.$fileName or die $!;
	foreach my $header (@head){
	    if($header =~/#PBS(.+)\.$number/){
		chomp $header;
		print NBPS $header.'.'.$counter."\n";
	    }else{
		print NBPS $header;
	    }
	}
	my $to = ($i+$step-1>$#commands?$#commands:$i+$step-1);
	my @slice = @commands[$i..$to];
	print NBPS @slice;
	print NBPS @tail;
	$counter++;
	close NBPS;
        #build a list of file names and then run all of them together
	push @files2sub, $fileName;
    }
}
foreach my $sub (@files2sub){
    print $sub,"\n";
    system "qsub $sub";
}
