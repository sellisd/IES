#!/usr/bin/perl
use warnings;
use strict;
# read multiple sequence alignment files
# sort based on number of sequencies
# create multiple .pbs files and qsub them
my $alnPath = $ARGV[0];
my $outPath = '/pandata/sellis/msas/trees/';
my $listfilesPath = '/pandata/sellis/msas/listfiles/';
opendir(DH, $alnPath) or die $!;
my @files = grep {/.*\.nucl\.fa/} readdir(DH);
close DH;
my @sizes;
foreach my $file (@files){
    my $filesize = -s $alnPath.$file;
    push @sizes, $filesize;
}
my @indices = sort {$sizes[$a] <=> $sizes[$b]} 0..$#sizes;
my $block = 50;
my $pbsCounter = 0;

my @filesSorted; # files sorted by increasing file size

foreach my $index (@indices){
    if ($sizes[$index]>0){ # exclude empty files
	push @filesSorted, $files[$index];
    }
}

for(my $i = 0; $i <= $#filesSorted; $i+=$block){
    # make PBS file
    open PBS, '>phyl.'.$pbsCounter.'.pbs' or die $!;
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N phyl.',$pbsCounter,"\n";
    print PBS '#PBS -e error.',$pbsCounter,"\n";
    print PBS '#PBS -o output.',$pbsCounter,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    open LST, '>'.$listfilesPath.'listfiles.'.$pbsCounter or die $!;
    for(my $j = 0; $j<$block; $j++){ #make listfiles
	my $file = $filesSorted[$i+$j];
	next unless defined $file;
	print LST $alnPath.$file,"\n";
    }
    close LST;
    #call Simon's script
    my $command = '/pandata/penel/exemple_raxml/ms_fasta_files_2_srs.pl  --option=\'-m GTRGAMMA\' --type=dna -modele=DNA -partition /pandata/sellis/msas/listfiles/listfiles.'.$pbsCounter.' '.$outPath.'out.'.$pbsCounter.'.srs';
    print PBS "$command","\n";
    close PBS;
    print 'qsub phyl.'.$pbsCounter.'.pbs'."\n";
    system 'qsub phyl.'.$pbsCounter.'.pbs';
    $pbsCounter++;
}
