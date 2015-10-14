#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Bio::TreeIO;
use Getopt::Long;

my $help;
my @ingroup;
my @outgroup;
my $S;
my $treePath;
my $output;
my $usage  = <<HERE;

Find subtrees rooted at a given speciation node and annotates the tips as as in- or out-groups

usage pickClade.pl [OPTIONS]
where OPTIONS can be:
-ingroup:  comma separated list of ingroup species (no space)
-outgroup: comma separated list of outgroup species names (no space)
-path:     path where PHYLDOG reconciled trees are located
-output:   output file name
-s:        
-help|?:   this help screen

e.g.

./pickClade.pl -s 0 -ingroup Paramecium_biaurelia,Paramecium_tetraurelia,Paramecium_sexaurelia -outgroup Paramecium_caudatum  -path ~/data/IES_data/msas/phyldog/results/ -output ~/data/IES_data/rdb/Pca.PbiPtePse.dat

./pickClade.pl -s 2 -ingroup Paramecium_biaurelia,Paramecium_tetraurelia -outgroup Paramecium_sexaurelia  -path ~/data/IES_data/msas/phyldog/results/ -output ~/data/IES_data/rdb/Pse.PbiPte.dat

HERE

    die $usage unless (GetOptions(
			   'help|?'     => \$help,
			   's=i'        => \$S,
			   'ingroup=s'  => \@ingroup,
			   'outgroup=s' => \@outgroup,
			   'path=s'     => \$treePath,
			   'output=s'   => \$output
		       ));

@ingroup = split(/,/, join(',', @ingroup));
@outgroup = split(/,/, join(',', @outgroup));
die $usage if $help;
die $usage if $#ingroup == -1 or $#outgroup == -1;
die $usage unless defined($S);

my %group;

foreach my $species (@ingroup){
    $group{$species} = 1;
}
foreach my $species (@outgroup){
    $group{$species} = -1;
}

opendir DH, $treePath or die $!;
my @files = grep {/^.*.\.ReconciledTree$/} readdir(DH);
close DH;

open OUT, '>'.$output or die $!;
print OUT "subtree\tcluster\tgeneId\tgroup\tspecies\n"; #header

my $subtreeId = 0; # counter of subtrees
my $fileCounter = 0;
foreach my $fileName(sort @files){
  $fileName =~ /(\d+)\.ReconciledTree$/;
  print $fileCounter, '/', $#files, "\r";
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
		    die "Multiple ingroups from the current speciation node. Perhaps the speciation node and the ingroup/outgroup were not defined properly";
		}else{
		    @ingroup = @leaves;
		}
	    }elsif ($group eq 'outgroup'){
		if(@outgroup){
		    die "Multiple outgroups from the current speciation node. Perhaps the speciation node and the ingroup/outgroup were not defined properly";
		}else{
		    @outgroup = @leaves;
		}
	    }else{
		die;
	    }
	}
	foreach my $geneId (@ingroup){
	    print OUT $subtreeId, "\t", $cluster, "\t", $geneId, "\t", 'ingroup', "\t", gene2species($geneId),"\n";
	}
	foreach my $geneId (@outgroup){
	    print OUT $subtreeId, "\t", $cluster, "\t", $geneId, "\t", 'outgroup', "\t", gene2species($geneId),"\n";
	}
	$subtreeId++;
    }
  }
  $fileCounter++;
}

close OUT;

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
	#  print $group, ' ', $element,"\n";
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
