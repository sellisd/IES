#!/usr/bin/perl
use warnings;
use strict;
use Bio::Tools::GFF;
use lib'.';
use functions;
use Getopt::Long;

my $usage = <<HERE;

Extract relevant IES information from a gff3 file and print them in a tab delimited file
usage iesInfo.pl file

HERE

my $gff = $ARGV[0];

my $gffI = Bio::Tools::GFF->new('-file' => $gff,
 				'-gff_version' => 3);
&printab('id', 'scaffold', 'altSeqNo', 'start', 'end', 'upstreamFlank', 'downstreamFlank', 'length', 'front', 'back', 'sequence');
while(my $feature = $gffI->next_feature()){
    my $scaffold = $feature->seq_id;
    my $start = $feature->start;
    my $end = $feature->end;
    my $id;
    my $name;
    my $sequence;
    my @altSequences;
    my $floating;
    my $species;
    my $age;
    my $upstreamFlank;
    my $downstreamFlank;
    my $length;
    my $lengthClass;
    my $front;
    my $back;
    my $addTA = 0;
    foreach my $tag ($feature->get_all_tags()){
	my @values = $feature->get_tag_values($tag);
	if ($tag ne 'alternative_seq' and $#values > 0){
	    die "$gff: more than one value per tag: $tag @values";
	}
	if ($tag eq 'Name'){
	    $name = $values[0];
	}elsif($tag eq 'ID'){
	    $id = $values[0];
	}elsif ($tag eq 'bounded_by_ta'){
	    $addTA = 1;
	}elsif ($tag eq 'sequence'){
	    $sequence = $values[0];
	}elsif ($tag eq 'junction_seq'){
	    $values[0] =~ /^([actgn]+)TA([actgn]+)$/ or die "$gff: error with regexp for junction sequence: $values[0]";
	    $upstreamFlank = uc($1).'TA';
	    $downstreamFlank = 'TA'.uc($2);
	    if(length($upstreamFlank) != length($downstreamFlank)){
		die "$gff: error with flank sequence sizes!";
	    }
	}elsif($tag eq 'alternative_seq'){
	    @altSequences = @values;
	}elsif($tag eq 'score'){
	    # ignore score
	}else{
	    die "$gff: unknown tag $tag";
	}
    }
    if($addTA == 1){
	$sequence .= 'TA';
	$length = length($sequence) - 2;
    }else{
	$length = length($sequence);
    }
    die "$gff name ($name) and id ($id) do not match" if $name ne $id;
    &printab($id, $scaffold, $#altSequences, $start, $end, $upstreamFlank, $downstreamFlank, $length, substr($sequence, 0, 20), substr($sequence, -20, 20), $sequence);
}
