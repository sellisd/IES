#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use Getopt::Long;

my $runNo;
my $toq;
my $path;
my $help;
my $step;

my $usage = <<HERE;

Check the results of an MSA run on the cluster and rerun those necessary by splitting independent
multiple sequence alignments and changing waiting queue.

usage: splitPBS.pl [OPTIONS]
where OPTIONS can be:
    -run:    number of run, 1 is the first (from msa.pl) and the rest increment
    -toq:    in which queue to submit the new jobs (q1hour, q1day or q1week)
    -path:   path where pbs files are located
    -step:   break each pbs in how many independent runs
    -help|?: this help screen

HERE

die $usage unless(GetOptions(
		      'help|?'   => \$help,
		      'toq=s'    => \$toq,
		      'run=i'    => \$runNo,
		      'step=i'   => \$step,
		      'path=s'   => \$path
		  ));

die $usage unless -d $path;
opendir(DH, $path) or die $!;
my @files;
if($runNo == 1){
    @files = grep{/error\.(\d+)*/} readdir(DH);
}elsif($runNo == 2){
    @files = grep{/error\.\d+\.(\d+)*/} readdir(DH);
}elsif($runNo == 3){
    @files = grep{/error\.\d+\.\d+\.(\d+)*/} readdir(DH);
}else{
    die "unexpected run number";
}
close DH;

my @files2sub;

foreach my $file (@files){
    next unless -s catfile($path, $file); #only error files with not zero size
    $file =~ /error.([.\d]+)/; # extract number
    my $number = $1;
    #if error file an expected warning
    open ERR, catfile($path, $file) or die $!;
    my $line = readline(ERR);
    chomp $line;
    if ($line =~ /^mkdir: cannot create directory(.*)File exists$/){
	next;
    }
    close ERR;
    open PBS, catfile($path,'msa.'.$number.'.pbs') or die $!;
    my @head;
    my @commands;
    my @tail;
    while(my $line = <PBS>){
	if (substr($line,0,1) eq '#'){
	    chomp $line;
	    if ($line =~ /#PBS\s+-q\s+(q1hour)|(q1day)|(q1week)/){
		$line = '#PBS -q '.$toq;
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
#    system "qsub $sub";
}
