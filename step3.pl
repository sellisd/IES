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
my $rbP = '/home/dsellis/tools/revbayes/projects/cmake/rb'; # revBayes path

my $nr = getNotation($notationF);
my $asrRevF = '/home/dsellis/projects/IES/src/asr.Rev';

# prepare revBayes runs
my $cmdl = 'Rscript --vanilla preRev.R > /home/dsellis/data/IES/analysis/preRev.log';
run($cmdl, 1);

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

    setOutAsr($asrRevF, $rbRun1Rev, '/pandata/sellis/asr'.$asrRun.'/', '/pandata/sellis/asr'.$asrRun.'/run1/', '/pandata/sellis/asr'.$asrRun.'/geneFamilies.dat', 1); # first time through print node index
    setOutAsr($asrRevF, $rbRun2Rev, '/pandata/sellis/asr'.$asrRun.'/', '/pandata/sellis/asr'.$asrRun.'/run2/', '/pandata/sellis/asr'.$asrRun.'/geneFamilies.dat', 0);
}

#move to the cluster

#    run("$rbP $rbRun1Rev", 1);
#    run("$rbP $rbRun2Rev", 1);

exit(1);
#=======================


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
