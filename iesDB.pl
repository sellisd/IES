#!/usr/bin/perl
use warnings;
use strict;
use Bio::Tools::GFF;
use File::HomeDir;
use File::Spec::Functions qw(catfile);
use lib'.';
use functions;
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
my $gff = '/home/dsellis/data/IES/analysis/filtscaf/poc.ies.gff3';
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
	    die "more than one value per tag: $tag @values";
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
	    $values[0] =~ /^([actgn]+)TA([actgn]+)$/ or die "error with regexp for junction sequence: $values[0]";
	    $upstreamFlank = uc($1).'TA';
	    $downstreamFlank = 'TA'.uc($2);
	    if(length($upstreamFlank) != length($downstreamFlank)){
		die "error with flank sequence sizes!";
	    }
	}elsif($tag eq 'alternative_seq'){
	    @altSequences = @values;
	}else{
	    die "unknown tag $tag";
	}
    }
    if($addTA == 1){
	$sequence .= 'TA';
	$length = length($sequence) - 2;
    }else{
	$length = length($sequence);
    }
    die if $name ne $id;
    &printab($id, $scaffold, $#altSequences, $start, $end, $upstreamFlank, $downstreamFlank, $length, substr($sequence, 0, 20), substr($sequence, -20, 20), $sequence);
}

# #extract raw information and create a "database" of IES interesting information
# #species id start end length sequence wgd?
# #INPUT ies.gff3 file
# #OUTPUT tab file
# #my $file = '/home/dsellis/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.fl.gff3';
# my $file = '/home/dsellis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3';
# print "id\tscaffold\tstart\tend\tlength\tsequence\twgd1\twgd2\twgd3\n"; # for P. tetraurelia
# open IN, $file or die $!;
# while(my $line = <IN>){
#     next if substr($line,0,1) eq '#';
#     chomp $line;
#     (my $scaffold, my $method, my $type, my $start, my $end,my $score, my $strand, my $phase, my $attributes) = split " ", $line;
#     my @tags = split ";", $attributes;
#     my $id;
#     my $name;
#     my $junction_seq;
#     my $sequence;
#     my $bounded_by_ta;
#     my $wgd1;
#     my $wgd2;
#     my $wgd3;
#     foreach my $tag (@tags){
# 	$tag =~ /(.*)=(.*)/;
# 	my $tagName = $1;
# 	my $tagValue = $2;
# 	if($tagName eq 'ID'){
# 	    $id = $tagValue;
# 	}elsif($tagName eq 'Name'){
# 	    $name = $tagValue;
# 	}elsif($tagName eq 'junction_seq'){
# 	    $junction_seq = $tagValue;
# 	}elsif($tagName eq 'sequence'){
# 	    $sequence = $tagValue;
# 	}elsif($tagName eq 'bounded_by_ta'){
# 	    $bounded_by_ta = $tagValue;
# 	}elsif($tagName eq 'conserved_IES_WGD1'){
# 	    $wgd1 = ($tagValue eq 'yes'?1:0)
# 	}elsif($tagName eq 'conserved_IES_WGD2'){
# 	    $wgd2 = ($tagValue eq 'yes'?1:0)
# 	}elsif($tagName eq 'conserved_IES_WGD3'){
# 	    $wgd3 = ($tagValue eq 'yes'?1:0)
# 	}
#     }
#     if ($bounded_by_ta){
# 	$sequence .= 'TA';
#     }
#     my $length = length($sequence);
# #    print "$id $start $end $length $sequence $junction_seq\n"; # for P. caudatum
#     print "$id\t$scaffold\t$start\t$end\t$length\t$sequence\t$wgd1\t$wgd2\t$wgd3\n"; # for P. tetraurelia
# }
# close IN;


# ies border db
# id length begin10bp end10bp
