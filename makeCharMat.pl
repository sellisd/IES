#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $help;
my $noAbsence = 0;

my $usage = <<HERE;

Prepare character matrices for multibiphy runs. The program reads the
filtered nucleotide alignments and .gblock files. The output depends
on the -noAbsence option. If FALSE prints the caracter matrices with
zeroes based on the .gblock files and create empty character
matrices. If TRUE prints only character matrices in phylipp format.
usage makeCharMat.pl [-noAbsence]

INTPUT
    /alignments/filtered/.aln                filtered protein alignments (>50% prot id)
    /alignments/filtered/.aln.fasta.gblocks  gblocks coordinates for the filtered alignments
OUTPUT
    /alignments/charMatphy/.phy             phylip format output character matrices

HERE

die $usage unless(GetOptions('help|?'    => \$help,
			     'noAbsence' => \$noAbsence
		  ));
die $usage if $#ARGV >1;
die $usage if $help;

my $filteredAlnPath = '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $gblocksPath =     '/home/dsellis/data/IES_data/msas/alignments/filtered/';
my $charMatPath =     '/home/dsellis/data/IES_data/msas/alignments/charMat/';
my $outCharMatPath =  '/home/dsellis/data/IES_data/msas/alignments/charMatphy/';

mkdir $outCharMatPath unless -d $outCharMatPath;
opendir(DH, $filteredAlnPath) or die $!; 
my @treeF = grep {/cluster\.(\d+)\.aln$/} readdir(DH); # find all clusters which passed the filtering

foreach my $fileName (@treeF){
    $fileName =~ /cluster\.(\d+)\.aln$/;
    my $cluster = $1;  # get the name of the cluster from the filename
    print $cluster,"\t";
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
	print "0\n";
	next; #skip files that have no conserved blocks
    }else{
	print "1\t";
    }
    # find if there is a character matrix (alignment has IES in gblocks)
    my $iesCharMatF = $charMatPath.'cluster.'.$cluster.'.dat';
    my $extCharMatF = $outCharMatPath.'cluster.'.$cluster.'.phy';
    my @matrix;
    my @genes;    
    my $lineCounter = 0;

    if(-e $iesCharMatF){ # if the cluster has IES
	print "1\n";
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
		push @ar, "0"x$extraColumns unless $noAbsence;
		my $geneName = shift @ar;
		my $string = join('',@ar);
		push @matrix, $geneName."\t".$string;
	    }
	    $lineCounter++;
	}
	close CM;
    }else{ # if not
	print "0\n";
	# from alignment find name of genes
	if($noAbsence){
	    next; #skip alignments with no IES
	}
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
#    print $cluster,"\n";use Data::Dumper;
#    print Dumper @matrix;
    open OUTM, '>'.$extCharMatF or die "$! $extCharMatF";
    print OUTM  $#matrix + 1,"\t",$alignmentLength,"\n";
    print OUTM join("\n",@matrix),"\n";
    close OUTM;
}
