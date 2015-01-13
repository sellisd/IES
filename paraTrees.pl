#!/usr/bin/perl
use warnings;
use strict;
# wrapper for trees.R, runs in parallel @pbil
my $inputPath = $ARGV[0];
opendir(DH, $inputPath) or die "$inputPath $!";
my @files = readdir(DH);
close DH;
foreach my $file (@files) {
    next unless -f $inputPath.$file;
    next unless $file =~ /.*\.nucl\.fa$/;
    my $group = $file;
    $group =~ s/cluster\.(\d+)\.nucl\.fa/$1/;
    my $pbsF = $group.'.pbs';
    open PBS, '>'.$pbsF or die $!;
    my $cmdString = 'Rscript --vanilla trees.R '.$inputPath.$file.' '.$group;
    #write pbs file
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N tree.',$group,"\n";
    print PBS '#PBS -e error.',$group,"\n";
    print PBS '#PBS -o output.',$group,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    print PBS "$cmdString","\n";
    close PBS;
    print "qsub ".$group.'.pbs',"\n";

    close PBS;
#make PBS file
#qsub
}
