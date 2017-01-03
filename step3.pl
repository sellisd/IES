#!/usr/bin/perl
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;

my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};

my $rbP = '/home/dsellis/tools/revbayes-1.0.0/projects/cmake/rb'; # revBayes path
my $notationF  = catfile($basePath, 'analysis', 'notation.csv');
my $nr         = getNotation($notationF);
my $logP       = catdir($basePath, 'analysis', 'log');
my $asrRevF    = './asr.Rev';

# prepare revBayes runs
my $cmdl = 'Rscript --vanilla preRev.R > '.$logP.'preRev.log';
run($cmdl, 1);

# test for low starting probability files

#./testRb.pl asr1 > testRb1.sh
#./testRb.pl asr2 > testRb2.sh
#./testRb.pl asr3 > testRb3.sh
#chmod 744 testRb?.sh
#./testRb1.sh
#./testRb2.sh
#./testRb3.sh

#3285   1  2  3
#5456   0  2  3
#5663   0  0  3
#10007  1  2  3
#11561  1  2  3

# manually remove low starting probability files from genefamilies

my %filtFam = (3285  => 1,
               5456  => 1,
               5663  => 1,
               10007 => 1,
               11561 => 1);
for(my $asrRun = 1; $asrRun <= 3; $asrRun++){
    my $gfFile = catfile($basePath, 'asr'.$asrRun, 'geneFamilies.dat');
    my @keep;
    open GF, $gfFile or die "$gfFile: $!";
    while (my $line = <GF>){
        chomp $line;
        unless($filtFam{$line}){
            push @keep, $line;
        }
    }
    close GF;
    open GF, '>', $gfFile or die "$gfFile: $!";
    print GF join("\n", @keep);
    close GF;
}

for(my $asrRun = 1; $asrRun<=3; $asrRun++){
    my $rbResultsP = catdir($basePath, 'analysis', 'asr'.$asrRun);
    # directories with rb output
    my $rbRun1         = catdir($rbResultsP, 'run1');
    my $rbRun2         = catdir($rbResultsP, 'run2');
    my $rbNodeIndexesP = catdir($rbResultsP, 'rbNodeIndexes');
    # Rev scripts
    my $rbRun1Rev = catfile($rbRun1, 'asr1.Rev');
    my $rbRun2Rev = catfile($rbRun2, 'asr2.Rev');
    make_path($rbResultsP)     unless -e $rbResultsP;
    make_path($rbNodeIndexesP) unless -e $rbNodeIndexesP;
    make_path($rbRun1)         unless -e $rbRun1;
    make_path($rbRun2)         unless -e $rbRun2;
    setOutAsr($asrRevF,
              $rbRun1Rev,
              catdir($basePath, 'asr'.$asrRun),
              catdir($basePath, 'asr'.$asrRun, 'run1'),
              catfile($basePath, 'asr'.$asrRun, 'geneFamilies.dat'),
              1); # first time through print node index
    setOutAsr($asrRevF,
              $rbRun2Rev,
              catdir($basePath, 'asr'.$asrRun),
              catdir($basePath, 'asr'.$asrRun, 'run2'),
              catfile($basePath, 'asr'.$asrRun, 'geneFamilies.dat'),
              0);

    run("$rbP $rbRun1Rev".' >'.$logP.'asr'.$asrRun.'run1.log'.' &', 0); #run both at the same time
    run("$rbP $rbRun2Rev".' >'.$logP.'asr'.$asrRun.'run2.log', 0);
}

sub setOutAsr{
    # modify the basic Rev script for multiple runs, and optionally keep specific lines
    my $file = shift @_;
    my $outfile = shift @_;

    my $indir = shift @_;
    my $outdir = shift @_;
    my $geneFamilyFile = shift @_;
    my $saveIndex = shift @_;
    open IN, $file or die $!;
    open OUT, '>', $outfile or die $!;
    while(my $line = <IN>){
	$line =~ s/INDIR/$indir/g;
	$line =~ s/OUTDIR/$outdir/g;
	$line =~ s/GENEFAMILYFILE/$geneFamilyFile/g;
	$line =~ s/SAVEINDEX/$saveIndex/g;
	print OUT $line;
    }
    close OUT;
    close IN;
}
