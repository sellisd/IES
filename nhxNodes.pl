#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;
use IO::String;
use Getopt::Long;

# read a nhx tree file, and for each node print the node ID and as a key the name of all leaves

my $help;
my $rb;
my $usage = <<HERE;

Reads nhx trees file from phyldog output and creates node key files. Files with the node name in the ND section and a string of the leaf names corresponding to that node concatenated in alphanumeric order. If the node is a leaf only its own name is included.
usage

nhxNodes.pl [OPTIONS] inputfiles

where options can be:
  rb      preprocess input from revBayes to make it [&&NHX]
  help|?  this help screen
HERE

die $usage unless(GetOptions('rb'     => \$rb,
                             'help|?' => \$help));
die $usage if $help;

foreach my $treeF (@ARGV){
    print $treeF,"\n";
    my $input;
    if($rb){   # if revBayes node index files preprocess and read from string
	open IN, $treeF or die $!;
	my $line = readline(IN);
	close IN;
	$line =~ s/\[&index=/[&&NHX:ND=/g;
	my $inputString = IO::String->new($line);
	$input = new Bio::TreeIO(-fh     => $inputString,
				 -format => "nhx");
    }else{
	$input = new Bio::TreeIO(-file   => $treeF,
				 -format => "nhx");
    }
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
