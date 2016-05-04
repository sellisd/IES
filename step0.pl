#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;

# directories
my $bamD = '/home/dsellis/data/IES/bam/';
my $bedD = '/home/dsellis/data/IES/analysis/bed';

#files
my $pcaCovFiltered = catfile($bedD, 'pca.cov.be'); # coverage over 14x
my $pcaCovFiltmerge = catfile($bedD, 'pca.covm.be'); # merged
# calculate regions of P. caudatum with genomic coverage > 14x
my $cmdl = 'bedtools genomecov -ibam '.
    catfile($bamD, 'PCAUD_MIC10_BCP_AFIOSF_2.BOWTIE.pc_mac_43c3d_v1.1.pe.sorted.bam').
    ' -bg | awk \'$4 > 14\' > '.
    $pcaCovFiltered;
run($cmdl, 0);

$cmdl = "bedtools merge -i $pcaCovFiltered > $pcaCovFiltmerge";
run($cmdl, 0);
