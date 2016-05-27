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
my $coverageD = '/home/dsellis/data/IES/analysis/coverage';
my $bedD      = '/home/dsellis/data/IES/analysis/bed';

make_path($coverageD) unless -d $coverageD;
my $nr = getNotation($notationF);

foreach my $sp (sort keys %$nr){
  my %pab = %{$nr->{$sp}}; #de-reference for less typing
  # calculate genome coverage
  my $coverageF = catfile($coverageD, $pab{'abr'}.'.cov');
  my $geneBedF  = catfile($bedD, $pab{'abr'}.'.gene.be');
  my $geneCovF = catfile($coverageD, $pab{'abr'}.'.gene.cov');
  my $cmdl = $bedtools.' genomecov '.
      ' -ibam '.catfile($pab{'datapath'}, $pab{'micbam'}).
      ' -g '.catfile($pab{'datapath'}, $pab{'MacF'}).
      ' -d > '.$coverageF;
  run($cmdl, 1);
  # and also genome coverage per scaffold

  # calculate average per nucleotide read coverage over genes
  run("./geneCov.pl -cov $coverageF -bed $geneBedF -out $geneCovF", 1);
}



# directories
# my $bamD = '/home/dsellis/data/IES/bam/';
# my $bedD = '/home/dsellis/data/IES/analysis/bed';

#files
# my $pcaCovFiltered = catfile($bedD, 'pca.cov.be'); # coverage
# my $genesF = catfile($bedD, 'pca.gene.be'); # genes

# calculate average per nulcotide coverage over genes
# my $cmdl = './geneCov.pl > /home/dsellis/data/IES/analysis/tables/pca.gene.cov';
# run($cmdl, 1);
# my $pcaCovFiltmerge = catfile($bedD, 'pca.covm.be'); # merged
# # calculate regions of P. caudatum with genomic coverage > 14x
# my $cmdl = 'bedtools genomecov -ibam '.
#     catfile($bamD, 'PCAUD_MIC10_BCP_AFIOSF_2.BOWTIE.pc_mac_43c3d_v1.1.pe.sorted.bam').
#     ' -bg | awk \'$4 > 14\' > '.
#     $pcaCovFiltered;
# run($cmdl, 0);

# $cmdl = "bedtools merge -i $pcaCovFiltered > $pcaCovFiltmerge";
# run($cmdl, 0);
