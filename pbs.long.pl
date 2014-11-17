#!/usr/bin/perl
use warnings;
use strict;
#change que to pbs runs that do not fit in 1 hour
my $path = $ARGV[0];
opendir(DH, $path) or die $!;
my @files = readdir(DH);
close DH;
foreach my $file (@files){
    next unless -s $path.$file;
    next unless $file =~ /error\.(\d+)\.(\d+)/; #keep only 2nd error files
    my $number = $1;
    my $secondNumber = $2;
    my $msaF = 'msa.'.$number.'.'.$secondNumber.'.pbs';
    open PBS, $msaF or die "$!: $msaF";
    my @slurp = <PBS>;
    my $contents = join '', @slurp;
    $contents =~ s/q1hour/q1day/m;
    $contents =~ s/(msa|error|output)\.($number)\.($secondNumber)/$1.$2.$3.long/mg;
    close PBS;
    my $newmsaF = $msaF;
    $newmsaF =~ s/\.pbs/.day.pbs/;
    open OUTPBS, '>'.$newmsaF or die $!;
    print OUTPBS $contents;
    close OUTPBS;
    my $cmdString = "qsub $newmsaF";
    print $cmdString,"\n";
    system $cmdString;
 }
