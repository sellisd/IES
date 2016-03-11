#!/usr/bin/perl
use warnings;
use strict;

opendir(DH, $ARGV[0]) or die $!;
my @files = grep{/.*\.pbs/} readdir(DH);
close DH;
my %errors;

foreach my $file (@files){
    #check if pbs file has a correctly formated output and an empty error file
    $file =~ /(.*).pbs/;
    my $chunk = $1;
    my $outputFile = 'output.'.$chunk;
    my $errorFile = 'error.'.$chunk;
    if(-e $outputFile){
	open OUT, $outputFile or die $!;
	my $lineCounter = 0;
	my $validOutput = 0;
	while(my $line = <OUT>){
	    chomp $line;
	    if ($lineCounter == 1){
		chomp $line;
		if ($line eq 'telos'){
		    $validOutput = 1;
		}
	    }
	    $lineCounter++;
	}
	close OUT;
	if ($validOutput == 0){
	    &addErrorMsg($chunk, 'Not valid output');
	}
    }else{
	&addErrorMsg($chunk, 'Not run');
    }
    if(-s $errorFile){
	open ER, $errorFile or die $!;
	my $lineCounter = 0;
	while (my $line = <ER>){
	    chomp $line;
	    if ($line =~/^Selenocysteine \(U\) at position \d+ replaced by X$/){
		# ignore warning
	    }elsif($line =~ /Warning: One or more U or O characters replaced by X for alignment score calculations at positions/){
		# ignore warning
	    }elsif($line eq ''){
		# ignore empty lines
	    }else{
		&addErrorMsg($chunk, 'Interrupted');
	    }
	    $lineCounter++;
	}
	close ER;
    }
}

# print summary of  errors

foreach my $chunks (keys %errors){
    print "$chunks, @{$errors{$chunks}}\n"
}

# prepare to rerun locally all the failed runs or rerun in another queue?

open RL, '>runlocal.sh' or die $!;
foreach my $run (keys %errors){
    print RL 'blastp -query /home/dsellis/data/IES/tempdat/fastachunks/chunk.'.$run.'.fa -db /home/dsellis/data/IES/analysis/protdb/allprot -outfmt 6 -out /home/dsellis/data/IES/analysis/allvsall/blastout.chunk.'.$run.'.fa.dat -seg yes > '.$run.'.log &',"\n";
}
close RL;

sub addErrorMsg{
    my $code = shift @_;
    my $message = shift @_;
    if(defined($errors{$code})){
	push @{$errors{$code}}, $message;
    }else{
	$errors{$code} = [$message];
    }
}
