#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;
use lib'.';
use functions;

# Homolog Finder finds groups of homologs (orthologs and paralogs) descending from a speciation node

my $speciationNode = 0;

my $treePath = '/home/dsellis/data/IES_data/msas/phyldog/results/';
my $matricesPath = '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';
opendir DH, $treePath or die $!;
my @files = grep {/^.*.\.ReconciledTree$/} readdir(DH);
close DH;

print "id\tcluster\tgeneId\tspecies\n"; #header
my $id = 0; # counter of homologous groups
foreach my $fileName(sort @files){
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    my $matrixF = $matricesPath.'cluster.'.$cluster.'.phy';
    if (-e $matrixF){
#	print $cluster,":\n";
	my $input = new Bio::TreeIO(-file   => $treePath.$fileName,
				    -format => "nhx");
	my $tree = $input->next_tree;
	foreach my $nodeO ($tree->get_nodes){
	    my $speciationEvent = $nodeO->get_tag_values('S');
	    my $eventType = $nodeO->get_tag_values('Ev');
	    if($eventType eq 'S' and $speciationEvent eq $speciationNode){
		my $found = 0;
		foreach my $descendentsO ($nodeO->get_all_Descendents){
		    if($descendentsO->is_Leaf()){
			$found = 1; # at least one leaf descends from the speciation node
			my $gene = $descendentsO->id();
			my $speciesName = &gene2species(substr($gene, 0, 4));
			print $id, "\t", $cluster, "\t", $gene, "\t", $speciesName, "\n";
		    }
		}
		if ($found == 1){
		    $id++;
		}else{
		    die "speciation node is a leaf";
		}
	    }
	}
    }else{
#	print $cluster," has no matrix\n";
    }
}
