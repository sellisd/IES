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
my $notationF =  catfile($homeD, 'data/IES/analysis/notation.tab');

open N, $notationF or die $!;
my $header = readline(N);
my %notation;
while(my $line = <N>){
    chomp $line;
(my $abbreviation, my $datapath, my $binomial, my $taxId, my $geneGff, my $cdsF, my $protF, my $geneF, my $MacF, my $iesGff, my $annotation, my $prefix) = split "\t", $line;
    $notation{$binomial} = {
	'annotation' => $annotation,
	'prefix'     => $prefix,
	'abr'        => $abbreviation,
	'datapath'   => $datapath,
	'taxId'      => $taxId,
	'geneGff'    => $geneGff,
	'cdsF'       => $cdsF,
	'protF'      => $protF,
	'geneF'      => $geneF,
	'MacF'       => $MacF,
        'iesGff'     => $iesGff
    };
}
close N;

# Variables
my $scafCutoff = 10_000;


# calculate scaffold lengths
# --------------------------
if(0){
system("./scaffoldStats.pl /home/dsellis/data/IES/primaurelia/genome/pprimaurelia_mac_AZ9-3_v1.0.fa        > /home/dsellis/data/IES/analysis/ppr.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/biaurelia/genome/biaurelia_V1-4_assembly_v1.fasta        > /home/dsellis/data/IES/analysis/pbi.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/tetraurelia/genome/ptetraurelia_mac_51.fa                > /home/dsellis/data/IES/analysis/pte.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/pentaurelia/genome/ppentaurelia_mac_87_v1.0.fa           > /home/dsellis/data/IES/analysis/ppe.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/sexaurelia/genome/sexaurelia_AZ8-4_assembly_v1.fasta     > /home/dsellis/data/IES/analysis/pse.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/octaurelia/genome/poctaurelia_mac_138_v1.0.fa            > /home/dsellis/data/IES/analysis/poc.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/tredecaurelia/genome/ptredecaurelia_209_AP38_filtered.fa > /home/dsellis/data/IES/analysis/ptr.scaf");
system("./scaffoldStats.pl /home/dsellis/data/IES/caudatum/caudatum_43c3d_assembly_v1.fasta >/home/dsellis/data/IES/analysis/pca.scaf");
}

# filter scaffolds
# ----------------

if(0){
    foreach my $sp (keys %notation){
	my $fsref = buildPaths($sp, $notation{$sp}{'annotation'});
	my @args = (
	    '-length',  catfile($homeD, $fsref->{'scaffoldF'}),
	    '-protein', catfile($homeD, $fsref->{'proteinF'}),
	    '-gff',     catfile($homeD, $fsref->{'anotF'}),
	    '-species', $fsref->{'pab'},
	    '-gene',    catfile($homeD, $fsref->{'geneF'}),
	    '-ies',     catfile($homeD, $fsref->{'iesF'}),
	    '-outdir',  catfile($homeD, '/data/IES/analysis/filtscaf/'),
	    '-cutoff',  $scafCutoff
	    );
	print "@args\n";
	system("./filterScaffolds.pl", @args);
    }
}

# make BLAST protein database
# ---------------------------

#find /home/dsellis/data/IES/analysis/ -name "*.protein.fa"
if(0){
    make_path("/home/dsellis/data/IES/analysis/protdb/");
system('cat /home/dsellis/data/IES/analysis/filtscaf/ptr.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/ppe.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/ppr.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/pca.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/poc.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/pse.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/pte.protein.fa \
     /home/dsellis/data/IES/analysis/filtscaf/pbi.protein.fa \
     /home/dsellis/data/IES/thermophila/gene/T_thermophila_June2014.protein.fa \
     >  /home/dsellis/data/IES/analysis/protdb/allprot.fa');
system('makeblastdb -in /home/dsellis/data/IES/analysis/protdb/allprot.fa -dbtype prot -out  /home/dsellis/data/IES/analysis/protdb/allprot');
}

# split fasta files for faster processing in parallel
# ---------------------------------------------------
if(0){
    system('./preblast.pl /home/dsellis/data/IES/analysis/protdb/allprot.fa /home/dsellis/data/IES/tempdat/fastachunks/');
}

# parse basic ies information
#----------------------------
if(0){
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/ppr.ies.gff3 > ~/data/IES/analysis/filtscaf/ppr.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/pbi.ies.gff3 > ~/data/IES/analysis/filtscaf/pbi.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/pte.ies.gff3 > ~/data/IES/analysis/filtscaf/pte.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/ppe.ies.gff3 > ~/data/IES/analysis/filtscaf/ppe.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/pse.ies.gff3 > ~/data/IES/analysis/filtscaf/pse.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/poc.ies.gff3 > ~/data/IES/analysis/filtscaf/poc.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/ptr.ies.gff3 > ~/data/IES/analysis/filtscaf/ptr.ies.tab";
    system "./iesInfo.pl ~/data/IES/analysis/filtscaf/pca.ies.gff3 > ~/data/IES/analysis/filtscaf/pca.ies.tab";
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
    system "Rscript --vanilla gene2intergenic.R ~/data/IES/analysis/bed/pca.gene.be ~/data/IES/analysis/bed/pca.inter.be ~/data/IES/analysis/pca.scaf";
}


# find overlap in bed files
if(0){
    foreach my $sp (keys %notation){
	my $abr = $notation{$sp}{'abr'};
	next if $abr eq 'pso';
	my $cmdl = 'bedtools intersect -a ~/data/IES/analysis/bed/'.$abr.'.ies.be'.
	    ' -b ~/data/IES/analysis/bed/'.$abr.'.cds.be'.
	    ' -b ~/data/IES/analysis/bed/'.$abr.'.intron.be'.
	    ' -b ~/data/IES/analysis/bed/'.$abr.'.inter.be'.
	    ' -names cds intron intergenic -wo > ~/data/IES/analysis/bed/'.$abr.'.IESin.be';
	print $cmdl, "\n";
	system $cmdl;
    }
}

# make genebank files and incorporate ies information

if(0){
    system "./makeGeneBank.pl -species ppr";
    system "./makeGeneBank.pl -species pbi";
    system "./makeGeneBank.pl -species pte";
    system "./makeGeneBank.pl -species ppe";
    system "./makeGeneBank.pl -species pse";
    system "./makeGeneBank.pl -species poc";
    system "./makeGeneBank.pl -species ptr";
    system "./makeGeneBank.pl -species pso";
    system "./makeGeneBank.pl -species pca";
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
