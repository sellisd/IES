#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;

#scan reconciled gene trees for pairs of orthologous IES between *P. tetraurelia* and *P. biaurelia* (their last common ancestor is a speciation event giving rise to *P. tetraurelia* and *P. biaurelia*).
# The number for this speciation event is hard coded to be 4

my $treePath = '/home/dsellis/data/IES_data/msas/phyldog/results/';
my $matricesPath = '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';
opendir DH, $treePath or die $!;
my @files = grep {/^.*.\.ReconciledTree$/} readdir(DH);
close DH;
foreach my $fileName(sort @files){
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    my $matrixF = $matricesPath.'cluster.'.$cluster.'.phy';
    if (-e $matrixF){
#	print $cluster,":\n";
	my $input = new Bio::TreeIO(-file   => $treePath.$fileName,
				    -format => "nhx");
	my $tree = $input->next_tree;
	my @leafObjects = $tree->get_leaf_nodes;
	my @bi;
	my @te;
	foreach my $leafO (@leafObjects){
	    my $gene = $leafO->id();
	    if (substr($gene, 0, 4) eq 'PBIA'){
		push @bi, $gene;
	    }
	    if(substr($gene,0,4) eq 'PTET'){
		push @te,$gene;
	    }
	}
	for(my $i = 0; $i <= $#bi; $i++){
	    my $pbiNode = $tree->find_node(-id => $bi[$i]);
	    for(my $j = $i; $j <= $#te; $j++){
		my $pteNode = $tree->find_node(-id => $te[$j]);
		my $lca = $tree->get_lca(-nodes => [$pbiNode, $pteNode]);
		my $speciationEvent = $lca->get_tag_values('S');
		my $eventType = $lca->get_tag_values('Ev');
		if($eventType eq 'S' and $speciationEvent == 4){
		    print $cluster, "\t", $speciationEvent, "\t", $bi[$i], "\t", $te[$i],"\n";
		}
	    }
	}  
    }else{
#       print $cluster," has no matrix\n";
    }
}
