#!/usr/bin/perl
use warnings;
use strict;
use Bio::TreeIO;
# parse PHYLDOG output remove leafs for which we have no IES information (P. caudatum and T. thermophila)
# Read character matrices, remove the corresponding lines and insert a number of all zero columns equal to the length of the aligned block

my $home = "/home/dsellis/";
my $treeDir = $home.'data/IES_data/msas/phyldog/result';
opendir DH, $treeDir or die $!;
my @files = grep {/.*\.ReconciledTree$/} readdir(DH);
close DH;
foreach my $fileName (@files){
    $fileName =~ /(\d+)\.ReconciledTree$/;
    my $cluster = $1;
    print $cluster,"\n";
    my $file = $home.'data/IES_data/msas/phyldog/result/'.$cluster.'.ReconciledTree';
    my $output = $home.'data/IES_data/msas/asr/trees/'.$cluster.'.tre';
    my $input = new Bio::TreeIO(-file => $file, -format => "nhx"); 
    my $tree = $input->next_tree;
    # find leaves to keep
    my %leavesToKeep;
    my @toKeepO;
    for my $node ($tree->get_leaf_nodes){
	my $id = $node->id;
	my $speciesAbr = substr($id,0,4);
	if($speciesAbr eq 'PCAU'
	   or $speciesAbr eq 'TTHE'){
	}else{
	    $leavesToKeep{$id} = 1;
	    push @toKeepO, $node;
	}
    }
    # move up the root to the MRCA of leaves with IES data
    my $newRoot = $tree->get_lca(\@toKeepO);
    $tree->set_root_node($newRoot);
    #    print $tree->as_text('newick'),"\n";die;
    # if there are remaining leaves without data remove them
    my @leavesToRemove;
    for my $node ($tree->get_leaf_nodes){
	my $id = $node->id;
	my $speciesAbr = substr($id,0,4);
	if($speciesAbr eq 'PCAU'
	   or $speciesAbr eq 'TTHE'){
	    push @leavesToRemove, $id;
	}
    }
    $tree->splice(-remove_id => \@leavesToRemove);
    my $out = new Bio::TreeIO(-file => ">$output",
			      -format => 'nhx');
    $out->write_tree($tree);

    # read gblocks file and count length of alignment

    my $gblockF = $home.'/data/IES_data/msas/alignments/aln/cluster.'.$cluster.'.nucl.fa.gblocks';
    my $alignmentLength = 0;
    open GB, $gblockF or die $!;
    while(my $line = <GB>){
	chomp $line;
	(my $cl, my $begin, my $end) = split " ", $line;
	$alignmentLength += $end - $begin + 1;
    }
    close GB;

    # read character matrix and add columns of zero and print to phylip format
    my $charMat = $home.'data/IES_data/msas/alignments/charMat/cluster.'.$cluster.'.L.dat';
    my $extCharMatF = $home.'data/IES_data/msas/asr/matrices/cluster.'.$cluster.'.phy';
    my @matrix;
    open OUTM, '>'.$extCharMatF or die $!;
    open CM, $charMat or die $!;
    my $lineCounter = 0;
    while(my $line = <CM>){
	chomp $line;
	my @ar = split " ", $line;
	my $extraColumns = $alignmentLength-$#ar;
	if (substr($ar[0],0,4) eq 'PCAU' or
	    substr($ar[0],0,4) eq 'TTHE'){
	    next; #skip rows
	}
	#print presence absence data
	if ($lineCounter == 0){
	    # do not print headers
	}else{
	    for(my $colCounter = 0; $colCounter <= $#ar; $colCounter++){
		if($colCounter > 0){
		    if ($ar[$colCounter] ne '0'){
			$ar[$colCounter] = '1';
		    }
		}
	    }
	    push @ar, "0"x$extraColumns;
	    my $geneName = shift @ar;
	    my $string = join('',@ar);
	    push @matrix, $geneName."\t".$string;
	}
	$lineCounter++;
    }
    print OUTM  $lineCounter-1,"\t",$alignmentLength,"\n";
    print OUTM join("\n",@matrix),"\n";
    close CM;
    close OUTM;
}

#parse phyldog and ancestral state reconstructions with ape (R)

# print the MPR trees
# parse species tree and keep names of nodes
#parse phyldog output
# MRCA
# 1. P2      P. biaurelia
# 2. P4      P. tetraurelia
# 3. P6      P. sexaurelia
# 4. Pc      P. caudatum
# 5. T       T. thermophila
# 6. P24     P. biaurelia & P. tetraurelia
# 7. P246    P. biaurelia & P. tetraurelia & P. sexaurelia
# 8. P246c   P. biaurelia & P. tetraurelia & P. sexaurelia & P. caudatum
# 9. P246cT  P. biaurelia & P. tetraurelia & P. sexaurelia & P. caudatum & T. thermophila

#The node numbers are hardcoded (data from OutputSpeciesTree_ConsensusNumbered.tree)

# my $MRCA1 = 8;  #  P. biaurelia
# my $MRCA2 = 7;  #  P. tetraurelia
# my $MRCA3 = 5;  #  P. sexaurelia
# my $MRCA6 = 6;  #  P. biaurelia & P. tetraurelia 
# my $MRCA7 = 4;  #  P. biaurelia & P. tetraurelia & P. sexaurelia
# my $asrF = $home.'data/IES_data/msas/alignments/asr/'.$cluster.'.dat';
#     my %hash;
#     open IN, $asrF or die "$! $asrF";
#     while (my $line = <IN>){
# 	chomp $line;
# 	next if $line eq 'lower upper';
# 	(my $node, my $lower, my $upper) = split " ", $line;
# 	$hash{$node} = [$lower,$upper];
#     }
#     close IN;
#for my $node ($tree->get_nodes){#
#	my @Ev = @{$node->{'_tags'}->{'Ev'}};
#	my @S = @{$node->{'_tags'}->{'S'}};
#	my @ND = @{$node->{'_tags'}->{'ND'}};
#	print $node->id,"\t";
#	print "Ev:@Ev S:@S ND:@ND\n";

#    print "   ",$node->bootstrap,"\n";
#}
#use Data::Dumper;
###print Dumper $tree;
