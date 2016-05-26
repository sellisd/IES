#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;
my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.csv');
my $rbResultsP = catfile($homeD, 'data/IES/analysis/asr/');

# directories with rb output
my $rbRun1 = catfile($rbResultsP, 'run1');
my $rbRun2 = catfile($rbResultsP, 'run2');
my $nodeIndexesP = catfile($rbResultsP, 'rbNodeIndexes');

# Rev scripts
my $rbRun1Rev = catfile($rbRun1, 'asr1.Rev');
my $rbRun2Rev = catfile($rbRun2, 'asr2.Rev');

my $nr = getNotation($notationF);
my $asrRevF = '/home/dsellis/projects/IES/src/asr.Rev';


my $cmdl = './charMats.pl > /home/dsellis/data/IES/analysis/iesdb/charMats.tab';
run($cmdl, 1);

$cmdl = 'Rscript --vanilla preRev.R > /home/dsellis/data/IES/analysis/preRev.log';
run($cmdl, 1);


make_path($rbRun1) unless -e $rbRun1;
make_path($rbRun2) unless -e $rbRun2;

setOutAsr($asrRevF, $rbRun1, $rbRun1Rev);
setOutAsr($asrRevF, $rbRun2, $rbRun2Rev);

run("rb $rbRun1Rev", 0, 1);
run("rb $rbRun2Rev", 1, 0);

run('./parseRevBayesAsr.pl > ~/data/IES/analysis/tables/avNodeProb.dat', 1);

# make a dictionary of node Ids across software
run('./nhxNodes.pl ~/data/IES/analysis/phyldog/results/*.ReconciledTree', 1);
run('./nhxNodes.pl -rb ~/data/IES/analysis/asr/rbNodeIndexes/nodeIndex.*.tre', 1);

sub setOutAsr{
    # modify the basic Rev script for multiple runs, and optionally keep specific lines
    my $file = shift @_;
    my $outpath = shift @_;
    my $outfile = shift @_;
    my $optline = shift @_; # should uncomment optional lines?
    open IN, $file or die $!;
    open OUT, '>', $outfile or die $!;
    while(my $line = <IN>){
	$line =~ s/OUTPATH/$outpath/g;
	if($optline){
	    $line =~ s/^#OPTLINE(.*)$/$1/;
	}
	print OUT $line;
    }
    close OUT;
    close IN;
}
