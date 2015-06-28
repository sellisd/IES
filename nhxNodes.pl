#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;

# read a nhx tree file, and for each node print the node ID and as a key the name of all leaves

my $treeF = $ARGV[0];

my $input = new Bio::TreeIO(-file   => $treeF,
                            -format => "nhx");
my $tree = $input->next_tree;


for my $node ( $tree->get_nodes ) {
    my @descNodes = $node->get_all_Descendents;
    my @leaves;
    foreach my $desc (@descNodes){
	my $id = $desc->{'_tags'}->{'ND'}[0];
	if($desc->is_Leaf){
	    my $name = $desc->id;
	    push @leaves,$name;
	}
	my $keyString;
	foreach my $leaf (sort @leaves){
	    $keyString .= '+'.$leaf;
	}
	print $id,"\t",$keyString,"\n";
    }

}
