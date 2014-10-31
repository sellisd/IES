#!/usr/bin
use warnings;
use strict;

opendir DH $ARGV[0] or die $!;
my @files = readir(DH);
close DH;

print `ls -1 output.* |wc`;
print `ls -1 error.* |wc`;
print `ls -1 *.pbs|wc`;
print `echo telos> 4comp; diff --from-file=4comp output.*; rm 4comp`;


foreach my $file (@files){
    next unless $file =~ /error.(\d+)/;
    open ER, $file or die $!;
    my $lineCounter = 0;
    while (my $line = <ER>){
        unless ($line =~ /^Selenocysteine (U) at position \d+ replaced by X$/){
            print "Error at $i line $lineCounter\n";
        }
        $lineCounter++;
    }
    close ER;
}
