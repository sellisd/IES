#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use File::Spec::Functions qw(catfile catdir);

my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};
my $tablesP  = catdir($basePath, 'analysis', 'tables');

# add homologous IES age to the iesdb by matching columns of the character matrix to the homologous IES id
my $asrRun = $ARGV[0];
my $linkF   = catfile($tablesP, 'homIES'.$asrRun.'.columns.link');
my $ageF    = catfile($tablesP, 'iesAge'.$asrRun.'.dat');
my $homiesF = catfile($tablesP, 'homIES.tab');

my %h; # gene family age information
my %l; # homologous IES id and column number link

# read name links
open LN, $linkF or die $!;
readline(LN); # header
while(my $line = <LN>){
    chomp $line;
    (my $geneFamily, my $id, my $column) = split " ", $line;
    $l{$geneFamily.'.'.$column} = $id;
}
close LN;
# read ages

open IN, $ageF or die $!;
readline(IN); # header
while(my $line = <IN>){
    chomp $line;
    (my $geneFamily, my $column, my $event) = split " ", $line;
    my $id = $l{$geneFamily.'.'.$column};
    die if defined($h{$geneFamily.'.'.$id});
    $h{$geneFamily.'.'.$id} = $event;
}
close IN;

# read database and fill in missing column
open DB, $homiesF or die "$! $homiesF";
my $header = readline(DB); # header
chomp $header;
print $header, "\t", 'age', "\n";
while(my $line = <DB>){
    chomp $line;
    (my $id, my $geneFamily, my $beginMSArange, my $endMSArange, my $gene, my $beginGene, my $endGene, my $beginMSA, my $endMSA, my $ies) = split " ", $line;
    if(defined($h{$geneFamily.'.'.$id})){
	print $line, "\t", $h{$geneFamily.'.'.$id}, "\n";
    }else{
	print $line, "\t", 'NA', "\n";
    }
}
close DB;
