#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $help;
my $burnIn = 100;
my $mbfPath = '/home/dsellis/data/IES_data/msas/asr/results/';

my $usage = <<HERE;

Read ancestral state reconstructions from multibiphy run and calculate summary probabilities ignoring the burnIn first steps.

parseMulitbiphy.pl [OPTIONS] outputFile

where OPTIONS can be:
 - burnIn:   Number of cycles to skip (default 100)
 - path:     Path to multibiphy treelist files (default : /home/dsellis/data/IES_data/msas/asr/results/)
 - help|?:   This help screen
HERE

die $usage unless (GetOptions('help|?'     => \$help,
			      'burnIn=i'   => \$burnIn,
			      'path=s'     => \$mbfPath
		   ));
die $usage if $help;
my $outputF;
if(defined($ARGV[0])){
    $outputF = $ARGV[0];
}else{
    print "No default output file\n";
    die $usage;
}

opendir DH, $mbfPath or die $!;
my @files = grep {/^ies\.cluster\.\d+\.phy\.treelist$/} readdir(DH);
close DH;

open OUT, '>', $outputF or die $!;
foreach my $fileName (@files){
#    $fileName = 'ies.cluster.10.phy.treelist';
    open IN, $mbfPath.$fileName or die $!;
    my $lineCounter = 0;
    my $totalCycles = 1;
    my @counts;
    my @prob;
    my $iesNo; #number of IES in cluster
    while(my $line = <IN>){
	chomp $line;
	my @ar = split " ", $line;
	if($lineCounter == 0){ #header
	    my @nodes = @ar;
	}else{
	    if($lineCounter > $burnIn){
		#in each column split by character
		my $colCounter = 0;
		foreach my $entry (@ar){
		    my @presenceAbsence = split "", $entry;
		    $iesNo = $#presenceAbsence + 1;
		    for(my $i = 0; $i <= $#presenceAbsence; $i++){
			if($presenceAbsence[$i] == 1){
			    if(defined($counts[$colCounter+$i])){
				$counts[2*$colCounter+$i]++;
			    }else{
				$counts[2*$colCounter+$i] = 1;
			    }
			}
		    }
		    $colCounter++;
		}
		$totalCycles++;
	    }
#array size is number of IES
	    #{node => [ies1Counter,ies2Counter,...]}
	}
	$lineCounter++;
    }
    close IN;
    $fileName =~ /^ies\.cluster\.(\d+)\.phy\.treelist$/;
    my $cluster = $1;
    print OUT $cluster,"\t$iesNo\t";
    foreach my $sum (@counts){
	if(defined($sum)){
	    push @prob, $sum/($totalCycles-1);
	}else{
	    push @prob, 0;
	}
    }
    print OUT join("\t",@prob),"\n";
}
close OUT;
# # for branch Si to Sj
# # calculate difference from all Si to Sj (including paralogs)
# # calculate difference from nodeA to nodeB (speciation nodes)
# my $nodeA = 0;
# my $nodeB = 1;
# # for all pair of paralogs compare Si to Sj (how similar they are)
# #read phyldog output find which nodes correspond to Si Sj of cluster X
# #read output of parseMultibify for cluster X calculate difference between node i and j for IES N
# my $phylTreePath = '/home/dsellis/data/IES_data/msas/phyldog1/results/';
# my $cluster = 10;
# open IN, $phylTreePath.$cluster.'.ReconciledTree' or die $!;
# my %spgH;  # hash with geneTreeNodes => species tree nodes
# my @fromNode;
# my @toNode;
# my $line = <IN>;
# while($line =~ /\[\&\&NHX:Ev=([SD]):S=(\d+):ND=(\d+)\]/g){
#     my $eventType = $1;
#     my $spNode    = $2;
#     my $node      = $3;
# #build correspondences (geneTreeNodes speciesNodes)
#     if($spNode == $nodeA){
# 	push @fromNode, $spNode;
#     }
#     if($spNode == $nodeB){
# 	push @toNode, $spNode;
#     }
#     if($eventType eq 'S'){
# 	$spgH{$node} = $spNode;
#     }
# }
# close IN;
# print "\n";
# foreach my $from (@fromNode){
#     print $prob[$from*$iesNo],"\t";
# }
# print "\n";
# foreach my $to (@toNode){
#     print $prob[$to*$iesNo],"\t";
# }
# print "\n";
# #chose comparisson , e.g 2->3 : Pbi,te -> Ptet
