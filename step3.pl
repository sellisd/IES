#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new(7);

my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.csv');
my $rbP = '/home/dsellis/tools/revbayes-1.0.0/projects/cmake/rb'; # revBayes path
my $logP = '/home/dsellis/data/IES/analysis/log/';
my $nr = getNotation($notationF);
my $asrRevF = '/home/dsellis/projects/IES/src/asr.Rev';

# base directory for output
#my $basedir = '/pandata/sellis/'; # on the cluster
my $basedir = '/home/dsellis/data/IES/analysis/'; # on the cluster


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
for(my $asrRun = 1; $asrRun<=3; $asrRun++){
    my $gfFile = catfile($basedir,'asr'.$asrRun.'/geneFamilies.dat') 
    open GF, $gfFile or die "$gfFile: $!"
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
    my $rbResultsP = catfile('/home/dsellis/data/IES/analysis/asr'.$asrRun.'/');
    # directories with rb output
    my $rbRun1 = catfile($rbResultsP, 'run1');
    my $rbRun2 = catfile($rbResultsP, 'run2');
    my $rbNodeIndexesP = catfile($rbResultsP, 'rbNodeIndexes');
    # Rev scripts
    my $rbRun1Rev = catfile($rbRun1, 'asr1.Rev');
    my $rbRun2Rev = catfile($rbRun2, 'asr2.Rev');
    make_path($rbResultsP) unless -e $rbResultsP;
    make_path($rbNodeIndexesP) unless -e $rbNodeIndexesP;
    make_path($rbRun1) unless -e $rbRun1;
    make_path($rbRun2) unless -e $rbRun2;
    setOutAsr($asrRevF, $rbRun1Rev, $basedir.'asr'.$asrRun.'/', $basedir.'asr'.$asrRun.'/run1/', $basedir.'asr'.$asrRun.'/geneFamilies.dat', 1); # first time through print node index
    setOutAsr($asrRevF, $rbRun2Rev, $basedir.'asr'.$asrRun.'/', $basedir.'asr'.$asrRun.'/run2/', $basedir.'asr'.$asrRun.'/geneFamilies.dat', 0);

    run("$rbP $rbRun1Rev".' >'.$logP.'asr'.$asrRun.'run1.log'.' &', 0); #run both at the same time
    run("$rbP $rbRun2Rev".' >'.$logP.'asr'.$asrRun.'run2.log', 0);
}

#move to the cluster

#    run("$rbP $rbRun1Rev", 1);
#    run("$rbP $rbRun2Rev", 1);




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
