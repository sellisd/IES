#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
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
    (my $binomial, my $annotation, my $prefix) = split "\t", $line;
    $notation{$binomial} = {
	'annotation' => $annotation,
	'prefix'     => $prefix
    };
}
close N;

# Variables
my $scafCutoff = 10_000;

sub buildPaths{
    my $binomial = shift @_;
    my $pabAnot = shift @_;
    $binomial =~ /^([A-Z])[a-z]+\s+([a-z]+)$/ or die;
    my $g = lc($1);
    my $spEpithet = $2;
    my $pab = $g.substr($spEpithet, 0, 2);
    my $scaffoldF = 'data/IES/analysis/'.$pab.'.scaf';
    my $proteinF  = 'data/IES/'.$spEpithet.'/gene/'.$pabAnot.'.protein.fa';
    my $geneF     = 'data/IES/'.$spEpithet.'/gene/'.$pabAnot.'.gene.fa';
    my $anotF     = 'data/IES/'.$spEpithet.'/gene/'.$pabAnot.'.gff3';
    my $iesF      = 'data/IES/'.$spEpithet.'/IES/'.$g.$spEpithet.'_internal_eliminated_sequence.gff3';
    return {
	'scaffoldF' => $scaffoldF,
	'proteinF'  => $proteinF,
	'geneF'     => $geneF,
	'anotF'     => $anotF,
	'iesF'      => $iesF,
	'pab'       => $pab
    }
}

# calculate scaffold lengths
# system("./scaffoldStats.pl /home/dsellis/data/IES/primaurelia/genome/pprimaurelia_mac_AZ9-3_v1.0.fa        > /home/dsellis/data/IES/analysis/ppr.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/biaurelia/genome/biaurelia_V1-4_assembly_v1.fasta        > /home/dsellis/data/IES/analysis/pbi.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/tetraurelia/genome/ptetraurelia_mac_51.fa                > /home/dsellis/data/IES/analysis/pte.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/pentaurelia/genome/ppentaurelia_mac_87_v1.0.fa           > /home/dsellis/data/IES/analysis/ppe.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/sexaurelia/genome/sexaurelia_AZ8-4_assembly_v1.fasta     > /home/dsellis/data/IES/analysis/pse.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/octaurelia/genome/poctaurelia_mac_138_v1.0.fa            > /home/dsellis/data/IES/analysis/poc.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/tredecaurelia/genome/ptredecaurelia_209_AP38_filtered.fa > /home/dsellis/data/IES/analysis/ptr.scaf");
# system("./scaffoldStats.pl /home/dsellis/data/IES/caudatum/caudatum_43c3d_assembly_v1.fasta >/home/dsellis/data/IES/analysis/pca.scaf");

# filter scaffolds

# foreach my $sp (keys %notation){
#     my $fsref = buildPaths($sp, $notation{$sp}{'annotation'});
#     my @args = (
# 	'-length',  catfile($homeD, $fsref->{'scaffoldF'}),
# 	'-protein', catfile($homeD, $fsref->{'proteinF'}),
# 	'-gff',     catfile($homeD, $fsref->{'anotF'}),
# 	'-species', $fsref->{'pab'},
# 	'-gene',    catfile($homeD, $fsref->{'geneF'}),
# 	'-ies',     catfile($homeD, $fsref->{'iesF'}),
# 	'-outdir',  catfile($homeD, '/data/IES/analysis/filtscaf/'),
# 	'-cutoff',  $scafCutoff
# 	);
#     print "@args\n";
#     system("./filterScaffolds.pl", @args);
# }


#find /home/dsellis/data/IES/analysis/ -name "*.protein.fa"
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
#`cat pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.protein.fa psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.protein.fa ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.protein.fa  pcaudatum_43c3d_annotation_v2.0/pcaudatum_43c3d_annotation_v2.0.protein.fa tthermophila/T_thermophila_June2014_proteins.fasta> working/combined.fa`


#`makeblastdb -in working/combined.fa -dbtype prot -out working/combined`
# with ncbi-blast+ version 2.2.30+
