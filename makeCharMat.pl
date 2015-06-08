#!/usr/bin/perl
use warnings;
use strict;
# read filtered nucleotide alignments
# read gblock files
# if present read character matrices and print to phylip format


my $filteredAlnPath = '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $gblocksPath =     '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $charMatPath =     '/home/dsellis/data/IES_data/msas/alignments/aln/charMat/';
my $outCharMatPath =  '';
opendir(DH, $filteredAlnPath) or die $!;
my @treeF = grep {/cluster\.(\d+)\.aln$/} readdir(DH);

foreach my $fileName (@treeF){
    $fileName =~ /cluster\.(\d+)\.aln$/;
    my $cluster = $1;
    print $cluster,"\n";
   # find Gblock blocks
    my $gblocksF = $gblocksPath.'cluster.'.$cluster.'.aln.fasta.gblocks';
    open IN, $gblocksF or die "$!: $gblocksF";
    my $alignmentLength = 0;
    while(my $line = <IN>){
	chomp $line;
	(my $cl, my $start, my $end) = split " ", $line;
	$alignmentLength += $end - $start + 1;
    }
    close IN;

    # find if there is a character matrix (alignment has IES in gblocks)
    my $iesCharMatF = $charMatPath.'cluster.'.$cluster.'.dat';
    my $extCharMatF = $charMatPath.'cluster.'.$cluster.'.phy';
    my @matrix;
    open OUTM, '>'.$extCharMatF or die "$! $extCharMatF";
    my @genes;    
    my $lineCounter = 0;
    if(-e $iesCharMatF){
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
    }else{
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
		push @matrix, $geneName."\t",$string;
		$lineCounter++;
	    }
	}
	close FASTA;
    }
    print substr($matrix[1],0,20),"\n";
    print OUTM  $lineCounter-1,"\t",$alignmentLength,"\n";
    print OUTM join("\n",@matrix),"\n";
    close OUTM;
}

