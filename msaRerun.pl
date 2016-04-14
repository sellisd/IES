#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use Getopt::Long;
use lib'.';
use functions;

my $in;
my $out;
my $toq;
my $path;
my $help;
my $step;

#first round: 
#./msaRerun.pl -toq q1hour -step 1 -in run1.dat -out run2.dat >msa2.log
#second round:
#./msaRerun.pl -toq q1day -step 1 -in run2.dat -out run3.dat >msa3.log
#third round:
#./msaRerun.pl -toq q1week -step 1 -in run3.dat -out run4.dat >msa4.log

my $usage = <<HERE;

Check the results of an MSA run on the cluster and rerun those necessary by splitting independent
multiple sequence alignments and changing waiting queue.

usage: msaRerun.pl [OPTIONS]
    where OPTIONS can be:
    -out:    file name to save failed run/geneFamily data
    -in:     file name with run/geneFamily data
    -toq:    in which queue to submit the new jobs (q1hour, q1day or q1week)
    -step:   break each pbs in how many independent runs
    -help|?: this help screen
HERE

die $usage unless(GetOptions(
		      'help|?'   => \$help,
		      'toq=s'    => \$toq,
		      'step=i'   => \$step,
		      'in=s'     => \$in,
		      'out=s'    => \$out
		  ));


my %runs;
open IN, $in or die $!;
$in =~ /run(\d+)\.dat/ or die $!;
my $runNo = $1;
open OUTT, '>run'.$runNo.'T.dat' or die $!;
while (my $line = <IN>){
    chomp $line;
    (my $run1, my $geneFamily) = split " ", $line;
    if(defined($runs{$run1})){
	next; # already tested
    }else{
	my $errorF = 'error.'.$run1;
	my $outputF = 'output.'.$run1;
	my $pbsF = 'msa.'.$run1.'.pbs';
	my $success = success({
	    'pbs'    => $pbsF,
	    'output' => $outputF,
	    'error'  => $errorF
			      });
	print OUTT $run1,"\t", $success,"\n";
	$runs{$run1} = $success;
    }
}
close OUTT;
close IN;

# If first round failed rerun in 2nd round
open my $outFH, '>', $out or die $!;
my @files2sub;
foreach my $didSucceed (sort keys %runs){
    if ($runs{$didSucceed}){
    }else{
	my $pbsF = 'msa.'.$didSucceed.'.pbs';
	my $f2s = splitPBS({
	    'pbs'   => $pbsF,
	    'toq'   => $toq,
	    'step'  => $step,
	    'outFH' => $outFH});
	push @files2sub, @$f2s;
    }

}
close $outFH;

foreach my $sub (@files2sub){
    print $sub,"\n";
    system "qsub $sub";
}

sub splitPBS{
    my $args  = shift @_;
    my $pbsF  = $args->{'pbs'};
    my $toq   = $args->{'toq'};
    my $step  = $args->{'step'};
    my $outFH = $args->{'outFH'};

    $pbsF =~ /msa\.([.\d]+)\.pbs/; # extract number
    my $number = $1;
    open PBS, $pbsF or die "$! $pbsF";
    my @head;
    my @commands;
    my @tail;
    my @files2sub;
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
 	my $fileName = 'msa.'.$number.'.'.$counter.'.pbs';
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
	foreach my $cmdl (@slice){
	    $cmdl =~ /cluster\.(\d+)\.fa/;
	    my $geneFamily = $1;
	    print $outFH $number.'.'.$counter,"\t", $geneFamily, "\n";
	}
 	print NBPS @slice;
 	print NBPS @tail;
 	$counter++;
 	close NBPS;
	# build a list of file names and then run all of them together
	push @files2sub, $fileName;
    }
    return \@files2sub;
}
