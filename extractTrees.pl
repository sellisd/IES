#!/usr/bin/perl
use warnings;
use strict;
#extract trees in newick format from RaxML output
my $path = $ARGV[0];
opendir(DH, $path) or die $!;
my @files = grep {/output\.*/} readdir(DH);
close DH;
foreach my $file (@files){
    open IN, $path.$file or die $!;
    my $cluster;
    while (my $line = <IN>){
	chomp $line;
	if (substr($line,0,7) eq 'Fichier'){
	    $line =~ /^.*cluster\.(\d+).nucl\.fa$/ or die $file;
	    $cluster = $1;
	}
	if (substr($line,0,6) eq 'Tree :'){
	    print 'cluster.'.$cluster,' ',substr($line,6),"\n";
	}
    }
}
close IN;
