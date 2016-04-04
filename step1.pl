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
    (my $abr, my $binomial, my $annotation, my $prefix) = split "\t", $line;
    $notation{$binomial} = {
	'annotation' => $annotation,
	'prefix'     => $prefix
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
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/ppr.ies.tab > ~/data/IES/analysis/filtscaf/ppr.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/pbi.ies.tab > ~/data/IES/analysis/filtscaf/pbi.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/pte.ies.tab > ~/data/IES/analysis/filtscaf/pte.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/ppe.ies.tab > ~/data/IES/analysis/filtscaf/ppe.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/pse.ies.tab > ~/data/IES/analysis/filtscaf/pse.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/poc.ies.tab > ~/data/IES/analysis/filtscaf/poc.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/ptr.ies.tab > ~/data/IES/analysis/filtscaf/ptr.ies.float";
    system "./annotateFloating.pl ~/data/IES/analysis/filtscaf/pca.ies.tab > ~/data/IES/analysis/filtscaf/pca.ies.float";
}


# check if we need to merge floating IES

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
