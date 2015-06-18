#!/usr/bin/perl
use warnings;
use strict;
# read filtered nucleotide alignments
# read gblock files
# if present read character matrices and print to phylip format
# if not create empty character matrix

# INPUT: /alignments/filtered/.aln                filtered protein alignments (>50% prot id)
#        /alignments/filtered/.aln.fasta.gblocks  gblocks coordinates for the filtered alignments

# OUTPUT: /alignments/charMatphy/.phy             phylip format output character matrices

my $filteredAlnPath = '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $gblocksPath =     '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $charMatPath =     '/home/dsellis/data/IES_data/msas/alignments/charMat/';
my $outCharMatPath =  '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';

opendir(DH, $filteredAlnPath) or die $!; 
my @treeF = grep {/cluster\.(\d+)\.aln$/} readdir(DH); # find all clusters which passed the filtering

foreach my $fileName (@treeF){
    $fileName =~ /cluster\.(\d+)\.aln$/;
    my $cluster = $1;  # get the name of the cluster from the filename
    print $cluster,"\n";
    # find Gblock blocks
    my $gblocksF = $gblocksPath.'cluster.'.$cluster.'.aln.fasta.gblocks';
    open IN, $gblocksF or die "$!: $gblocksF";
    my $alignmentLength = 0;
    while(my $line = <IN>){  # find how many columns are in the gblock
	chomp $line;
	(my $cl, my $start, my $end) = split " ", $line;
	$alignmentLength += $end - $start + 1;
    }
    close IN;
    if($alignmentLength == 0){
	next; #skip files that have no conserved blocks
    }
    # find if there is a character matrix (alignment has IES in gblocks)
    my $iesCharMatF = $charMatPath.'cluster.'.$cluster.'.dat';
    my $extCharMatF = $outCharMatPath.'cluster.'.$cluster.'.phy';
    my @matrix;
    open OUTM, '>'.$extCharMatF or die "$! $extCharMatF";

    my @genes;    
    my $lineCounter = 0;

    if(-e $iesCharMatF){ # if the cluster has IES
	open CM, $iesCharMatF or die $!;
	while(my $line = <CM>){
	    chomp $line;
	    my @ar = split " ", $line;
	    my $extraColumns = $alignmentLength-$#ar;
	    if (substr($ar[0],0,4) eq 'PCAU' or
		substr($ar[0],0,4) eq 'TTHE'){
		next; #skip rows
	    }
	    if($lineCounter == 0){
		# ignore header
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
	close CM;
    }else{ # if not
	# from alignment find name of genes
	open FASTA, $filteredAlnPath.'cluster.'.$cluster.'.nucl.fa' or die $!;
	while(my $line = <FASTA>){
	    if(substr($line,0,1) eq '>'){
		chomp $line;
		$line =~ /^>(.*)$/;
		my $geneName = $1;
		if (substr($geneName,0,4) eq 'PCAU' or
		    substr($geneName,0,4) eq 'TTHE'){
		    next; #skip rows
		}
		my $string = join('',"0"x$alignmentLength);
		push @matrix, $geneName."\t".$string;
		$lineCounter++;
	    }
	}
	close FASTA;
	if ($#matrix==-1){ # if there are clusters (without IES) that have only P. caudatum or T. thermophila
	    next;
	}
    }
    print OUTM  $#matrix + 1,"\t",$alignmentLength,"\n";
    print OUTM join("\n",@matrix),"\n";
    close OUTM;
}

