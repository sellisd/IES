#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;

# prepare character matrices
my $homIESdb = '/home/dsellis/data/IES/analysis/iesdb/homIESdb.dat';
open IN, $homIESdb or die $!;
readline(IN); # header
my %charMats;
while(my $line = <IN>){
    chomp $line;
    (my $id, my $geneFamily, my $beginMSArange, my $endMSArange, my $gene, my $beginGene, my $endGene, my $beginMSA, my $endMSA, my $ies) = split " ", $line;
    if(defined($charMats{$geneFamily}{$id}{$gene})){
	push @{$charMats{$geneFamily}{$id}{$gene}{'ies'}}, $ies;
	# data validation
	die unless ($charMats{$geneFamily}{$id}{$gene}{'begin'} == $beginMSArange);
	die unless ($charMats{$geneFamily}{$id}{$gene}{'end'} == $endMSArange);
    }else{
	$charMats{$geneFamily}{$id}{$gene} = {'begin' => $beginMSArange,
					      'end'   => $endMSArange,
					      'ies'   => [$ies]};
    }
}
close IN;
#use Data::Dumper; print Dumper %charMats;
# character matrix output
printab(('cluster', 'column', 'geneId', 'begin', 'end', 'ies'));
foreach my $geneFamily (sort {$a<=>$b} keys %charMats){
    foreach my $id (sort {$a <=> $b} keys %{$charMats{$geneFamily}}){
	foreach my $gene (sort keys  %{$charMats{$geneFamily}{$id}}){
	    printab($geneFamily, $id, $gene, $charMats{$geneFamily}{$id}{$gene}{'begin'},
		    $charMats{$geneFamily}{$id}{$gene}{'end'},
		    join(',', @{$charMats{$geneFamily}{$id}{$gene}{'ies'}}));
	}
    }
}
