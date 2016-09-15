#!/usr/bin/perl
use warnings;
use strict;
use File::Path qw(make_path);
use Bio::SeqIO;
use lib'.';
use functions;
use Getopt::Long;


my $inputFiltered = '/home/dsellis/data/IES/analysis/msas/filtered/';
my $outputPath = '/home/dsellis/data/IES/analysis/phyldog/';
##### dbg
use List::Util qw(shuffle);
$outputPath = '/home/dsellis/data/IES/analysis/phyldogPilot/';
make_path($outputPath) unless -d $outputPath;

#### dbg
my $fastaPath = 'aln/';
my $linkPath = 'link/';

#my $treePath = 'TreeFiles/';

#my $speciesTree = '(Tetrahymena_thermophila:5,(Paramecium_caudatum:0.3,(Paramecium_sexaurelia:0.15,(Paramecium_tetraurelia:0.1,Paramecium_biaurelia:0.1)):0.2)):0.01;';
#my $speciesTree = '(Tetrahymena_thermophila,(Paramecium_caudatum,(Paramecium_sexaurelia,(Paramecium_tetraurelia,Paramecium_biaurelia))));';

#make required folders
make_path($outputPath.$fastaPath) unless -d $outputPath.$fastaPath;
make_path($outputPath.$linkPath) unless -d $outputPath.$linkPath;
make_path($outputPath.'results') unless -d $outputPath.'results';
make_path($outputPath.'run') unless -d $outputPath.'run';


# open ST, '>'.$outputPath.'speciesTree.tre' or die;
# print ST $speciesTree."\n";
# close ST;

#read all trees in memory
#my %treesH;
#open TRE, $inputTreesF or die $!;
#while(my $line = <TRE>){
#    chomp $line;
#    (my $treeName, my $treeString) = split " ", $line;
#    my $cluster = substr($treeName,8,);
#    $treesH{$cluster} = $treeString;
#}
#close TRE;

#use filtered clusters
opendir(DH,$inputFiltered) or die $!;

my @charMats = grep {/cluster\.\d+\.nucl\.fa/} readdir(DH);
### dbg
@charMats = (shuffle(@charMats))[0..250-1];
###
foreach my $file (@charMats){
    $file =~ /cluster\.(\d+)\.nucl\.fa/;
    my $cluster = $1;
    my $fastaFileSource = $inputFiltered.'cluster.'.$cluster.'.nucl.fa';
    my $IF = Bio::SeqIO->new('-file'   => $fastaFileSource,
			     '-format' => 'Fasta');
    my %linkH;
    my @seqObjects;
    my %seqS;
    while(my $seqO = $IF->next_seq()){
	my $geneName = $seqO->display_id();
	push @seqObjects, $seqO;
	$seqS{$seqO->seq()} = 1;
	my $speciesName  = &gene2species($geneName);
	if(defined($linkH{$speciesName})){
	    push @{$linkH{$speciesName}}, $geneName;
	}else{
	    $linkH{$speciesName} = [$geneName];
	}
#	$OF->write_seq($seqO);
    }
    if(keys %seqS == 1){
	print "skipping $cluster\n";
	next;
    }
    my $fastaFileTarget = $outputPath.$fastaPath.$cluster.'.fasta'; 
    my $linkFile = $outputPath.$linkPath.$cluster.'.link';
    my $OF = Bio::SeqIO->new('-file'   => '>'.$fastaFileTarget,
			     '-format' => 'Fasta');
    foreach my $seqO (@seqObjects){
	$OF->write_seq($seqO);
    }
    open LIN, '>'.$linkFile or die $!;
    foreach my $speciesName (sort keys %linkH){
	print LIN $speciesName,':';
	print LIN join(';', @{$linkH{$speciesName}}),"\n";	
    }
    close LIN;
}
close DH;


#prepare appropriate directories and files for a phyldog run

#make 3 directories
# FastaFiles
#   cluster.X.fasta
# LinkFiles
#   cluster.X.link
#      speciesName:geneName1;geneName2;
# TreeFiles
#
