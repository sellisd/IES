#!/usr/bin/perl
use warnings;
use strict;
use Bio::Tools::GFF;
use lib'.';
use functions;
use Getopt::Long;

my $usage = <<HERE;

Extract information for genes from gff3 and prints .be files for CDS, introns and itergenic regions
usage iesInfo.pl file

usage: geneInfo.pl gffFile outBaseName
HERE



my $gff = $ARGV[0];
my $outbn = $ARGV[1];
my $cdsBEF = $outbn.'.cds.be';
my $geneBEF = $outbn.'.gene.be';
my $exonBEF = $outbn.'.exon.be';
open cdsF, '>'.$cdsBEF or die $!;
open geneF, '>'.$geneBEF or die $!;
open exonF, '>'.$exonBEF or die $!;

my $gffI = Bio::Tools::GFF->new('-file' => $gff,
 				'-gff_version' => 3);

while(my $feature = $gffI->next_feature()){
    my $scaffold = $feature->seq_id;
    my $start = $feature->start;
    my $end = $feature->end;
    my $type = $feature->primary_tag();
    my $id;
    my $name;
    my $sequence;
    my $parent;
    foreach my $tag ($feature->get_all_tags()){
	my @values = $feature->get_tag_values($tag);
	if ($tag eq 'Name'){
	    $name = $values[0];
	}elsif($tag eq 'ID'){
	    $id = $values[0];
	}elsif($tag eq 'Parent'){
	    $parent = $values[0];
	}else{
#	    die "unknown tag: $tag";
	}
    }
    die "$gff name ($name) and id ($id) do not match" if $name ne $id;
    if ($type eq 'CDS'){
	print cdsF "$scaffold\t$start\t$end\t$name\n";
    }elsif($type eq 'gene'){ # save gene begin end to get the intergenic complement
	print geneF "$scaffold\t$start\t$end\t$name\n";
    }elsif($type eq 'exon'){ # save the exon to get the complement introns
	print exonF "$scaffold\t$start\t$end\t$name\t$parent\n";
    }
}
    
close cdsF;
close geneF;
close exonF;
