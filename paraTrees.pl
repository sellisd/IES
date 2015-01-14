#!/usr/bin/perl
use warnings;
use strict;
# wrapper for trees.R
my $inputPath = $ARGV[0];
my $outputPath = $ARGV[1];
my $cluster = 0; # TRUE if run on the cluster
opendir(DH, $inputPath) or die "$inputPath $!";
my @files = readdir(DH);
close DH;
foreach my $file (@files) {
    next unless -f $inputPath.$file;
    next unless $file =~ /.*\.nucl\.fa$/;
    my $group = $file;
    $group =~ s/cluster\.(\d+)\.nucl\.fa/$1/;
    my $output = $file;
    $output =~ s/\.nucl.fa/.tre/;
    my $pbsF = $group.'.pbs';
    my $cmdString = 'Rscript --vanilla trees.R '.$inputPath.$file.' '.$group.' '.$outputPath.$output;
    if($cluster == 1){
	open PBS, '>'.$pbsF or die $!;
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
    }else{
	print $cmdString,"\n";
	system $cmdString;
    }
#make PBS file
#qsub
}
