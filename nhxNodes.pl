#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;
use Getopt::Long;

# read a nhx tree file, and for each node print the node ID and as a key the name of all leaves

my $help;
my $usage = <<HERE;

Reads nhx trees file from phyldog output and creates node key files. Files with the node name in the ND section and a string of the leaf names corresponding to that node concatenated in alphanumeric order. If the node is a leaf only its own name is included.
usage

nhxNodes.pl inputfiles

HERE

die $usage unless(GetOptions('help|?' => \$help));
die $usage if $help;

foreach my $treeF (@ARGV){
    print $treeF,"\n";
    my $input = new Bio::TreeIO(-file   => $treeF,
				-format => "nhx");
    my $tree = $input->next_tree;
    my $outF = $treeF.'.key';
    open OUT, '>'.$outF or die $!;
    for my $node ( $tree->get_nodes ) {
	my $keyString;
	my $id = $node->get_tag_values('ND');
	if($node->is_Leaf){ # if node is leaf
	    $keyString = $node->id;
	}else{
	    my @descNodes = $node->get_all_Descendents;
	    my @leaves;
	    foreach my $desc (@descNodes){
		if($desc->is_Leaf){
		    push @leaves,$desc->id;
		}
	    }
	    foreach my $leaf (sort @leaves){
		$keyString .= $leaf;
	    }
	}
    print OUT $id,"\t", $keyString,"\n";
    }
    close OUT;
}
