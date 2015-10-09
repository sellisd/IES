#!/usr/bin/perl
use warnings;
use strict;
use File::Path qw(make_path);
use Bio::SeqIO;
use lib'.';
use functions;
my $inputFiltered = '/home/dsellis/data/IES_data/msas/alignments/filtered/'; 
#my $inputTreesF = '/home/dsellis/data/IES_data/msas/tree/trees.tre';

my $outputPath = '/home/dsellis/data/IES_data/msas/phyldog/';
my $fastaPath = 'aln/';
my $linkPath = 'link/';

#my $treePath = 'TreeFiles/';

# do not use branch lengths
#my $speciesTree = '(Tetrahymena_thermophila:5,(Paramecium_caudatum:0.3,(Paramecium_sexaurelia:0.15,(Paramecium_tetraurelia:0.1,Paramecium_biaurelia:0.1)):0.2)):0.01;';
my $speciesTree = '(Tetrahymena_thermophila,(Paramecium_caudatum,(Paramecium_sexaurelia,(Paramecium_tetraurelia,Paramecium_biaurelia))));';

#make required folders
make_path($outputPath.$fastaPath) unless -d $outputPath.$fastaPath;
make_path($outputPath.$linkPath) unless -d $outputPath.$linkPath;
make_path($outputPath.'results') unless -d $outputPath.'results';
make_path($outputPath.'run') unless -d $outputPath.'run';
#make_path($outputPath.$treePath) unless -d $outputPath.$treePath;

open ST, '>'.$outputPath.'speciesTree.tre' or die;
print ST $speciesTree."\n";
close ST;

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
foreach my $file (@charMats){
    $file =~ /cluster\.(\d+)\.nucl\.fa/;
    my $cluster = $1;
    print $cluster,"\n";
    my $fastaFileSource = $inputFiltered.'cluster.'.$cluster.'.nucl.fa';
    my $fastaFileTarget = $outputPath.$fastaPath.$cluster.'.fasta'; 
    my $linkFile = $outputPath.$linkPath.$cluster.'.link';
#    my $treeFile = $outputPath.$treePath.$cluster.'.tre';
    my $IF = Bio::SeqIO->new('-file'   => $fastaFileSource,
			     '-format' => 'Fasta');
    my $OF = Bio::SeqIO->new('-file'   => '>'.$fastaFileTarget,
			     '-format' => 'Fasta');
    my %linkH;
    while(my $seqO = $IF->next_seq()){
	my $geneName = $seqO->display_id();
	my $speciesName  = &gene2species($geneName);
	if($speciesName eq 'Tetrahymena_thermophila'){
	    next;
	}
	if(defined($linkH{$speciesName})){
	    push @{$linkH{$speciesName}}, $geneName;
	}else{
	    $linkH{$speciesName} = [$geneName];
	}
	$OF->write_seq($seqO);
    }
#foreach sequence
# filter name
    # open FAI, $fastaFileSource or die $!;
    # open FAO, '>'.$fastaFileTarget or die $!;
    # while(my $line = <FAI>){
    # 	print FAO $line;
    # 	chomp $line;
    # 	if (substr($line,0,1) eq '>'){
    # 	    my $geneName = substr($line,1);
    # 	    my $speciesName = &gene2species($geneName);
    # 	    if(defined($linkH{$speciesName})){
    # 		push @{$linkH{$speciesName}}, $geneName;
    # 	    }else{
    # 		$linkH{$speciesName} = [$geneName];
    # 	    }
    # 	}
    # }
    # close FAI;
    # close FAO;
    open LIN, '>'.$linkFile or die $!;
    foreach my $speciesName (sort keys %linkH){
	print LIN $speciesName,':';
	print LIN join(';', @{$linkH{$speciesName}}),"\n";	
    }
    close LIN;
#    open TROUT, '>'.$treeFile or die $!; # <')(((((<
#    if(defined($treesH{$cluster})){
#	print TROUT $treesH{$cluster};
#    }else{
#        # In the case of a 3 gene tree the topology is a star
#	my @leafs;
#	foreach my $geneName(sort keys %linkH){
#	    push @leafs, @{$linkH{$geneName}};
#	}
#	print TROUT '(',join(',',@leafs),');',"\n";
#    }
#    close TROUT;
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
