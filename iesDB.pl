#!/usr/bin/perl
use warnings;
use strict;

#extract raw information and create a "database" of IES interesting information
#species id start end length sequence wgd?
#INPUT ies.gff3 file
#OUTPUT tab file
#my $file = '/home/dsellis/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.fl.gff3';
my $file = '/home/dsellis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3';
print "id\tscaffold\tstart\tend\tlength\tsequence\twgd1\twgd2\twgd3\n"; # for P. tetraurelia
open IN, $file or die $!;
while(my $line = <IN>){
    next if substr($line,0,1) eq '#';
    chomp $line;
    (my $scaffold, my $method, my $type, my $start, my $end,my $score, my $strand, my $phase, my $attributes) = split " ", $line;
    my @tags = split ";", $attributes;
    my $id;
    my $name;
    my $junction_seq;
    my $sequence;
    my $bounded_by_ta;
    my $wgd1;
    my $wgd2;
    my $wgd3;
    foreach my $tag (@tags){
	$tag =~ /(.*)=(.*)/;
	my $tagName = $1;
	my $tagValue = $2;
	if($tagName eq 'ID'){
	    $id = $tagValue;
	}elsif($tagName eq 'Name'){
	    $name = $tagValue;
	}elsif($tagName eq 'junction_seq'){
	    $junction_seq = $tagValue;
	}elsif($tagName eq 'sequence'){
	    $sequence = $tagValue;
	}elsif($tagName eq 'bounded_by_ta'){
	    $bounded_by_ta = $tagValue;
	}elsif($tagName eq 'conserved_IES_WGD1'){
	    $wgd1 = ($tagValue eq 'yes'?1:0)
	}elsif($tagName eq 'conserved_IES_WGD2'){
	    $wgd2 = ($tagValue eq 'yes'?1:0)
	}elsif($tagName eq 'conserved_IES_WGD3'){
	    $wgd3 = ($tagValue eq 'yes'?1:0)
	}
    }
    if ($bounded_by_ta){
	$sequence .= 'TA';
    }
    my $length = length($sequence);
#    print "$id $start $end $length $sequence $junction_seq\n"; # for P. caudatum
    print "$id\t$scaffold\t$start\t$end\t$length\t$sequence\t$wgd1\t$wgd2\t$wgd3\n"; # for P. tetraurelia
}
close IN;
