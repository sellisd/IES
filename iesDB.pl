#!/usr/bin/perl
use warnings;
use strict;

#extract raw information and create a "database" of IES interesting information
#species id start end length sequence wgd?
#INPUT ies.gff3 file
#OUTPUT tab file
my $file = '/home/dsellis/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.fl.gff3';
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
	}
    }
    if ($bounded_by_ta){
	$sequence .= 'TA';
    }
    my $lenght = length($sequence);
    print "$id $start $end $lenght $sequence $junction_seq\n";
}
close IN;
