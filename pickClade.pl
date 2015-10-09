#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Bio::TreeIO;
#given a speciation node number print in a reconciled tree the leafs of the node
# given an output group G and species A1-An of an input group pick subclades of a tree that have a monophyletic clade with the ingroup species A1-An and an outgroup G
# as PHYLDOG has already calculated the nodes corresponding to speciation events we can pick all the subtrees of this speciation node

my $treePath = '/home/dsellis/data/IES_data/msas/phyldog/results/';
my $matricesPath = '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';
my @ingroup = qw/Paramecium_biaurelia Paramecium_tetraurelia Paramecium_sexaurelia/;
my @outgroup = qw/Paramecium_caudatum/;
my %group;
my $S = 0; #number of speciation event
foreach my $species (@ingroup){
    $group{$species} = 1;
}
foreach my $species (@outgroup){
    $group{$species} = -1;
}

opendir DH, $treePath or die $!;
my @files = grep {/^.*.\.ReconciledTree$/} readdir(DH);
close DH;

print "cluster\tgeneId\tgroup\tspecies\n"; #header
foreach my $fileName(sort @files){
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    my $input = new Bio::TreeIO(-file   => $treePath.$fileName,
				-format => "nhx");
    my $tree = $input->next_tree;
    for my $nodeO ( $tree->get_nodes ) {
	my $speciationEvent = $nodeO->get_tag_values('S');
	my $eventType = $nodeO->get_tag_values('Ev');
	if($eventType eq 'S' and $speciationEvent == $S){
	    my @outgroup;
	    my @ingroup;
	    for my $child ($nodeO->each_Descendent){
		my @leaves;
		# should either be ingroup or outgroup
		if($child->is_Leaf){
		    push @leaves, scalar($child->id);
		}else{
		    for my $subTreeNode ($child->get_all_Descendents){
			if ($subTreeNode->is_Leaf){
			    push @leaves, scalar($subTreeNode->id);
			}
		    }
		}
		my $group =  &assignGroup(\@leaves);
		if ($group eq 'ingroup'){
		    if(@ingroup){
			die;
		    }else{
			@ingroup = @leaves;
		    }
		}elsif ($group eq 'outgroup'){
		    if(@outgroup){
			die;
		    }else{
			@outgroup = @leaves;
		    }
		}else{
		    die;
		}
	    }
	    foreach my $geneId (@ingroup){
		print $cluster, "\t", $geneId, "\t", 'ingroup', "\t", gene2species($geneId),"\n";
	    }
	    foreach my $geneId (@outgroup){
		print $cluster, "\t", $geneId, "\t", 'outgroup', "\t", gene2species($geneId),"\n";
	    }
	}
    }
}


#tidy output should be:
# cluster geneId group    species
# 100     PGM..  ingroup  PTET
# 100     MIC..  outgroup PCAU
sub assignGroup{
    my $arref = shift @_;
    my $counter = 0;
    my $group = 0; # positive ingroup negative outgroup
    foreach my $element (@$arref){
	$group += &isIngroup($element);
	$counter++;
#	print $group, ' ', $element,"\n";
    }
    if($counter == $group){
	return 'ingroup';
    }elsif($counter == -$group){
	return 'outgroup';
    }else{
	die;
    }
}

sub isIngroup{
    my $name = shift @_;
    my $species = gene2species($name);
    if(defined($group{$species})){
	return $group{$species};
    }else{
	die;
    }
}
