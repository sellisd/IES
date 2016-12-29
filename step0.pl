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

# binaries
my $bedtools = '/usr/local/bin/bedtools';

# output directories
my $coverageD = catfile($homeD, 'data/IES/analysis/coverage');
my $bedD      = catfile($homeD, 'data/IES/analysis/bed');
my $figuresD  = catfile($homeD, 'data/IES/analysis/figures/');

make_path($coverageD) unless -d $coverageD;
my $nr = getNotation($notationF);

foreach my $sp (sort keys %$nr){
  my %pab = %{$nr->{$sp}}; #de-reference for less typing
  # calculate genome coverage
  my $coverageF = catfile($coverageD, $pab{'abr'}.'.cov');
  my $avcovF = catfile($coverageD, $pab{'abr'}.'.avcov');
  my $geneBedF  = catfile($bedD, $pab{'abr'}.'.gene.be');
  my $geneCovF = catfile($coverageD, $pab{'abr'}.'.gene.cov');
  my $cmdl = $bedtools.' genomecov '.
      ' -ibam '.catfile($pab{'datapath'}, $pab{'micbam'}).
      ' -g '.catfile($pab{'datapath'}, $pab{'MacF'}).
      ' -d > '.$coverageF;
  run($cmdl, 0);
  # calculate average per nucleotide read coverage over genes
  run("./geneCov.pl -cov $coverageF -bed $geneBedF -out $geneCovF -avcov $avcovF", 0);
}



make_path($figuresD) unless -d $figuresD;
run("Rscript --vanilla ./geneFilterByCoverage.R", 0);