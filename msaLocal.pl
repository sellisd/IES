#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use File::Spec::Functions qw(catdir catfile);
use Parallel::ForkManager;
use File::Path qw(make_path);
use Time::HiRes qw(gettimeofday tv_interval);

# use mafft --auto for performing locally all the alignments
my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};

my $fastaPath = catdir($basePath, 'analysis', 'msas', 'fasta');
my $outPath   = catdir($basePath, 'analysis', 'msas', 'aln');
my $logPath   = catdir($basePath, 'analysis', 'msas', 'log');

make_path($outPath) unless -d $outPath;
make_path($logPath) unless -d $logPath;
opendir DH, $fastaPath or die "$! $fastaPath";
my @files = grep {/^cluster\.\d+\.fa$/} readdir(DH);

my $totalTime = 0;
foreach my $file (@files){
    my $outFile = $file;
    my $logFile = $file;
    $outFile =~ s/^(cluster\.\d+)\.fa$/$1.aln.fa/ or die $!;
    $logFile =~ s/^(cluster\.\d+)\.fa$/$1.log/ or die $!;
    my $inFP = catfile($fastaPath, $file);
    my $outFP = catfile($outPath, $outFile);
    my $logFP = catfile($logPath, $logFile);
    my $now = [gettimeofday];
    my $cmdl = "mafft --anysymbol --auto $inFP 1> $outFP 2> $logFP";
    $file =~/^cluster\.(\d+)\.fa$/;
    my $geneFamily = $1;
    print $geneFamily, ' ';
    system $cmdl;
    my $elapsed = tv_interval($now);
    print $elapsed, "\n";
    $totalTime += $elapsed;
}
close DH;

print "total: $totalTime sec for $#files gene families\n";
