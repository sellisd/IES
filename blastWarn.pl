#!/usr/bin/perl
use warnings;
use strict;

opendir(DH, $ARGV[0]) or die $!;
my @files = readdir(DH);
close DH;

print `ls -1 output.* |wc`;
print `ls -1 error.* |wc`;
print `ls -1 *.pbs|wc`;
#print `echo telos> 4comp; diff --from-file=4comp output.*; rm 4comp`;
my %output;
my %pbs;

foreach my $file (@files){
    if ($file =~ /output\.(\d+)/){
	$output{$1} = $file;
	open OUT, $file or die $!;
	my $lineCounter = 0;
	while(my $line = <OUT>){
	    if ($lineCounter == 1){
		chomp $line;
		print "Warning, strange output at $file\n" unless $line = 'telos';
	    }
	    $lineCounter++;
	}
	close OUT;
    }elsif($file =~ /(\d+).pbs/){
	$pbs{$1} = $file;
    }
    next unless $file =~ /error\.(\d+)/;
    open ER, $file or die $!;
    my $lineCounter = 0;
    while (my $line = <ER>){
	chomp $line;
	if ($line =~/^Selenocysteine \(U\) at position \d+ replaced by X$/){
	}elsif($line =~ /Warning: One or more U or O characters replaced by X for alignment score calculations at positions/){
	}elsif($line eq ''){
	}else{
	    print "Error at $file line $lineCounter\n";
        }
        $lineCounter++;
    }
    close ER;
}
foreach my $run (keys %pbs){
    unless(defined($output{$run})){
	print "Warning: ", $run, " did not finish\n";
    }
}
