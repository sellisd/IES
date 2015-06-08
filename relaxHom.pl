#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(all);

my $path = '/home/dsellis/data/IES_data/msas/alignments/aln/';
my $pathOut = '/home/dsellis/data/IES_data/msas/alignments/charMat/';

opendir(DH, $path) or die $!;
my @charMF = grep {/\.F\.dat$/} readdir(DH);
mkdir $pathOut unless -d $pathOut;

# read 2merge files and character matrices
# read gblocks and filter
foreach my $file (@charMF){
    $file =~ /cluster\.(\d+)/;
    my $cluster = $1;
    my @M;
    my @rowNames;
    my @colNames;
    open CM, $path.$file or die $!;
    my $rowCounter = 0;
    while(my $line = <CM>){
	chomp $line;
	my @ar = split " ", $line;
	if($rowCounter == 0){
	    @colNames = @ar;
	}else{
	    push @rowNames, shift @ar;
	    push @M, \@ar;
	}
	$rowCounter++;
    }
    close CM;
    my @NM;
    my $mergeF = 'cluster.'.$cluster.'.2merge.bed';
    my @NMcolNames;
    if(-s $path.$mergeF){
    #read 2merge file and build NM table
	open IN, $path.$mergeF or die $!;
	my $colCor = 0; #keep track of the NM and M column correspondence
	while(my $line = <IN>){
	    chomp $line;
	    (my $cl, my $start, my $end, my $cols) = split " ", $line;
	    my @merged = split ',', $cols;
	    push @NMcolNames, $start.'.'.$end;
	    if ($#merged > 0){
		for(my $i = 0; $i <= $#M; $i++){
		    my @mergedString;
		    foreach my $col2M (@merged){
			my $entry = $M[$i][$col2M];
			unless($entry eq '0'){
			    push @mergedString, $entry;
			}
		    }
		    if($#mergedString == -1){ #if all entries were 0
			push @mergedString, 0;
		    }
		    my $mergedStringS;
 		    if(all {$_ eq 'NA'} @mergedString){
			$mergedStringS = 'NA';
		    }else{
			$mergedStringS = join (',',@mergedString);
		    }
                    #new matrix gets merged column entries in place of the first of the merged ones
		    $NM[$i][$merged[0]-$colCor] = $mergedStringS;
		}
		$colCor+=$#merged;
	    }else{
		for(my $i = 0; $i <= $#M; $i++){
		    $NM[$i][$cols-$colCor] = $M[$i][$cols];#new matrix has the same column
		}
	    }
	}
	close IN;
	open OUT, '>'.$pathOut.'cluster.'.$cluster.'.dat' or die $!;
	print "Cluster $cluster\n";
	print OUT "geneName\t",join("\t",@NMcolNames),"\n";
	for(my $i = 0; $i <= $#NM; $i++){ # new matrix has the same number of rows
	    print OUT $rowNames[$i],"\t";
	    for(my $j = 0; $j <= $#{$NM[0]}; $j++){
		if(!defined($NM[$i][$j])){ #columns that our outside of Gblocks are not defined
		}else{
		    print OUT $NM[$i][$j],"\t";
		}
	    }
	    print OUT "\n";
	}
	close OUT;
    }
}
