#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;

# Evolution of IES in *Paramecium*

my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.csv');

my $nr = getNotation($notationF);

my %notation;

# Variables
my $scafCutoff = 10_000;

my $tablesP = catfile($homeD, 'data/IES/analysis/tables');

# calculate scaffold lengths
# --------------------------
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    my $genome = catfile($pab{'datapath'}, $pab{'MacF'});
    my $scafF = catfile($pab{'datapath'}, 'analysis/filtscaf', $pab{'abr'}.'.scaf');
    my $cmdl = './scaffoldStats.pl '.$genome.' > '.$scafF;
    run($cmdl, 1);
}

# filter scaffolds
# ----------------

foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $scafF   = catfile('/home/dsellis/data/IES/analysis/filtscaf', $pab{'abr'}.'.scaf');
    my $protein = catfile($pab{'datapath'}, $pab{'protF'}),
    my $gff     = catfile($pab{'datapath'}, $pab{'geneGff'}),
    my $gene    = catfile($pab{'datapath'}, $pab{'geneF'}),
    my $ies     = catfile($pab{'datapath'}, $pab{'iesGff'}),
    my $outdir  = catfile($pab{'datapath'}, 'analysis/filtscaf/'),
    my $cutoff  = $scafCutoff;
    my $cmdl = './filterScaffolds.pl'.
	' -length '.$scafF.
	' -ies '.$ies.
	' -protein '.$protein.
	' -gene '.$gene.
	' -gff '.$gff.
	' -cutoff '.$cutoff.
	' -species '.$pab{'abr'}.
	' -outdir '.$outdir;
    run($cmdl, 1);
}


# make BLAST protein database
# ---------------------------

#find /home/dsellis/data/IES/analysis/ -name "*.protein.fa"

my $protdbP = '/home/dsellis/data/IES/analysis/protdb/';
my $allprotF = catfile($protdbP, 'allprot.fa');
my $tthGenome = '/home/dsellis/data/IES/thermophila/gene/tth.protein.fa';
make_path($protdbP) unless -d $protdbP;
my @proteinF;
for my $sp (keys %$nr){
    my %pab = %{$nr->{$sp}};
    push @proteinF, catfile($pab{'datapath'}, 'analysis','filtscaf',$pab{'abr'}.'.protein.fa');
}
push @proteinF, $tthGenome; #add T. thermophila
my $cmdl = 'cat '.join(' ', @proteinF).' > '.$allprotF;
run($cmdl, 1);
$cmdl = 'makeblastdb -in '.$allprotF.' -dbtype prot -out '.catfile($protdbP,'allprot');
run($cmdl, 1);

# split fasta files for faster processing in parallel
# ---------------------------------------------------
run('./preblast.pl /home/dsellis/data/IES/analysis/protdb/allprot.fa /home/dsellis/data/IES/tempdat/fastachunks/', 1);

# parse basic ies information
#----------------------------
for my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $cmdl = './iesInfo.pl ~/data/IES/analysis/filtscaf/'.$pab{'abr'}.'.ies.gff3 > ~/data/IES/analysis/filtscaf/'.$pab{'abr'}.'.ies.tab';
    run($cmdl, 1);
}

# find which IES are floating
if(0){
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/ppr.ies.tab -floatF ~/data/IES/analysis/filtscaf/ppr.ies.float -bedF ~/data/IES/analysis/bed/ppr.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/pbi.ies.tab -floatF ~/data/IES/analysis/filtscaf/pbi.ies.float -bedF ~/data/IES/analysis/bed/pbi.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/pte.ies.tab -floatF ~/data/IES/analysis/filtscaf/pte.ies.float -bedF ~/data/IES/analysis/bed/pte.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/ppe.ies.tab -floatF ~/data/IES/analysis/filtscaf/ppe.ies.float -bedF ~/data/IES/analysis/bed/ppe.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/pse.ies.tab -floatF ~/data/IES/analysis/filtscaf/pse.ies.float -bedF ~/data/IES/analysis/bed/pse.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/poc.ies.tab -floatF ~/data/IES/analysis/filtscaf/poc.ies.float -bedF ~/data/IES/analysis/bed/poc.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/ptr.ies.tab -floatF ~/data/IES/analysis/filtscaf/ptr.ies.float -bedF ~/data/IES/analysis/bed/ptr.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/pso.ies.tab -floatF ~/data/IES/analysis/filtscaf/pso.ies.float -bedF ~/data/IES/analysis/bed/pso.ies.be";
    system "./annotateFloating.pl -tabF ~/data/IES/analysis/filtscaf/pca.ies.tab -floatF ~/data/IES/analysis/filtscaf/pca.ies.float -bedF ~/data/IES/analysis/bed/pca.ies.be";
}

# check if we need to merge floating IES
if(0){
    system "mkdir ~/data/IES/analysis/bed/";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/ppr.gff ~/data/IES/analysis/bed/ppr &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/pbi.gff ~/data/IES/analysis/bed/pbi &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/pte.gff ~/data/IES/analysis/bed/pte &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/ppe.gff ~/data/IES/analysis/bed/ppe &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/pse.gff ~/data/IES/analysis/bed/pse &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/poc.gff ~/data/IES/analysis/bed/poc &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/ptr.gff ~/data/IES/analysis/bed/ptr &";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/pso.gff ~/data/IES/analysis/bed/pso";
    system "./geneInfo.pl ~/data/IES/analysis/filtscaf/pca.gff ~/data/IES/analysis/bed/pca";
}

# create bed files
if(0){
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/ppr.exon.be ~/data/IES/analysis/bed/ppr.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/pbi.exon.be ~/data/IES/analysis/bed/pbi.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/pte.exon.be ~/data/IES/analysis/bed/pte.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/ppe.exon.be ~/data/IES/analysis/bed/ppe.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/pse.exon.be ~/data/IES/analysis/bed/pse.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/poc.exon.be ~/data/IES/analysis/bed/poc.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/ptr.exon.be ~/data/IES/analysis/bed/ptr.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/pso.exon.be ~/data/IES/analysis/bed/pso.intron.be";
    system "Rscript --vanilla exon2intron.R ~/data/IES/analysis/bed/pca.exon.be ~/data/IES/analysis/bed/pca.intron.be";
}

if(0){
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/ppr.gene.be ~/data/IES/analysis/bed/ppr.inter.be ~/data/IES/analysis/ppr.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/pbi.gene.be ~/data/IES/analysis/bed/pbi.inter.be ~/data/IES/analysis/pbi.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/pte.gene.be ~/data/IES/analysis/bed/pte.inter.be ~/data/IES/analysis/pte.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/ppe.gene.be ~/data/IES/analysis/bed/ppe.inter.be ~/data/IES/analysis/ppe.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/pse.gene.be ~/data/IES/analysis/bed/pse.inter.be ~/data/IES/analysis/pse.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/poc.gene.be ~/data/IES/analysis/bed/poc.inter.be ~/data/IES/analysis/poc.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/ptr.gene.be ~/data/IES/analysis/bed/ptr.inter.be ~/data/IES/analysis/ptr.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/pso.gene.be ~/data/IES/analysis/bed/pso.inter.be ~/data/IES/analysis/pso.scaf";
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/pca.gene.be ~/data/IES/analysis/bed/pca.inter.be ~/data/IES/analysis/pca.scaf";
}

# find overlap in bed files
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $cmdl = 'bedtools intersect -a ~/data/IES/analysis/bed/'.$abr.'.ies.be'.
	' -b ~/data/IES/analysis/bed/'.$abr.'.cds.be'.
	' -b ~/data/IES/analysis/bed/'.$abr.'.intron.be'.
	' -b ~/data/IES/analysis/bed/'.$abr.'.inter.be'.
	' -names cds intron intergenic -wo > ~/data/IES/analysis/bed/'.$abr.'.IESin.be';
    run($cmdl, 1);
}

# make IES table for iesDB
my $mergeF = '~/data/IES/analysis/filtscaf/ies2merge.dat';
$cmdl = "./reannotateFloating.pl ~/data/IES/analysis/filtscaf/*.ies.float > $mergeF";
run($cmdl, 1);

foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $cmdl = './iesdbTable.pl -merge '.$mergeF.
	' -iesin ~/data/IES/analysis/bed/'.$abr.'.IESin.be'.
	' -iestab ~/data/IES/analysis/filtscaf/'.$abr.'.ies.tab'.
	' -float ~/data/IES/analysis/filtscaf/'.$abr.'.ies.float'.
	' > ~/data/IES/analysis/iesdb/'.$abr.'.iesdb';
    run($cmdl, 1);
}

# make CDS table for iesDB
my $silixout = '/home/dsellis/data/IES/analysis/allvsall/blastout/silix.output';
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $cmdl = './cdsdbTable.pl -silixout '.$silixout.
	' -gff '.catfile($pab{'datapath'}, $pab{'geneGff'}).
	' > ~/data/IES/analysis/iesdb/'.$abr.'.cdsdb';
    run($cmdl, 1);
}

# find transcript coordinates for IES in genes
make_path($tablesP) unless -d $tablesP;
foreach my $sp (keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $cmdl = 'Rscript --vanilla ./iesInGenes.R '.
	'~/data/IES/analysis/bed/'.$abr.'.IESin.be '.
	'~/data/IES/analysis/iesdb/'.$abr.'.cdsdb '.
	catfile($tablesP, $abr.'.iesInGenes');
    run($cmdl, 1);
}
die;
exit(0);
#./maleTables.pl 
# make genebank files and incorporate ies information

die;
for my $sp (keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $cmdl = './makeGeneBank.pl -species '.$pab{'abr'};
    run($cmdl, 1);
}

if(0){
    system "./makeBed.pl ~/data/IES/analysis/gnbk/ppr.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/pbi.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/pte.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/ppe.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/pse.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/poc.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/ptr.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/pso.IES.gnbk";
    system "./makeBed.pl ~/data/IES/analysis/gnbk/pca.IES.gnbk";
}

if(0){
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/ppr.CDS.bed -b ~/data/IES/analysis/gnbk/ppr.IES.bed -wo > ~/data/IES/analysis/gnbk/ppr.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/pbi.CDS.bed -b ~/data/IES/analysis/gnbk/pbi.IES.bed -wo > ~/data/IES/analysis/gnbk/pbi.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/pte.CDS.bed -b ~/data/IES/analysis/gnbk/pte.IES.bed -wo > ~/data/IES/analysis/gnbk/pte.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/ppe.CDS.bed -b ~/data/IES/analysis/gnbk/ppe.IES.bed -wo > ~/data/IES/analysis/gnbk/ppe.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/pse.CDS.bed -b ~/data/IES/analysis/gnbk/pse.IES.bed -wo > ~/data/IES/analysis/gnbk/pse.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/poc.CDS.bed -b ~/data/IES/analysis/gnbk/poc.IES.bed -wo > ~/data/IES/analysis/gnbk/poc.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/ptr.CDS.bed -b ~/data/IES/analysis/gnbk/ptr.IES.bed -wo > ~/data/IES/analysis/gnbk/ptr.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/pso.CDS.bed -b ~/data/IES/analysis/gnbk/pso.IES.bed -wo > ~/data/IES/analysis/gnbk/pso.overlap";
    system "bedtools intersect -a ~/data/IES/analysis/gnbk/pca.CDS.bed -b ~/data/IES/analysis/gnbk/pca.IES.bed -wo > ~/data/IES/analysis/gnbk/pca.overlap";
}
