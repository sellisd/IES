#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use Parallel::ForkManager;
use File::Path qw(make_path);
use Time::HiRes qw(gettimeofday tv_interval);

# use mafft --auto for performing locally all the alignments
my $pm = Parallel::ForkManager->new(1);

my $fastaPath = '/home/dsellis/data/IES/analysis/msas/fasta/';
my $outPath = '/home/dsellis/data/IES/analysis/msas/aln/';
my $logPath = '/home/dsellis/data/IES/analysis/msas/log/';

make_path($outPath) unless -d $outPath;
make_path($logPath) unless -d $logPath;
opendir DH, $fastaPath or die "$! $fastaPath";
my @files = grep {/^cluster\.\d+\.fa$/} readdir(DH);

my @randomSample;
for(my $i = 0; $i <= 150; $i++){
    push @randomSample, $files[rand @files];
}
my $totalTime = 0;
foreach my $file (@randomSample){
#    my $pid = $pm->start and next;
    my $outFile = $file;
    my $logFile = $file;
    $outFile =~ s/^(cluster\.\d+)\.fa$/$1.aln.fa/ or die $!;
    $logFile =~ s/^(cluster\.\d+)\.fa$/$1.log/ or die $!;
    my $inFP = catfile($fastaPath, $file);
    my $outFP = catfile($outPath, $outFile);
    my $logFP = catfile($logPath, $logFile);
    my $now = [gettimeofday];
    my $cmdl = "mafft --anysymbol --auto $inFP 1> $outFP 2> $logFP";
    system $cmdl;
    my $elapsed = tv_interval($now);
    print $elapsed, "\n";
    $totalTime += $elapsed;

 #   $pm->finish;
}
close DH;

print "total: $totalTime sec for $#randomSample gene families\n";
