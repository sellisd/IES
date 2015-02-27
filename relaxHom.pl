#!/usr/bin/perl
use warnings;
use strict;

my $path = '/home/dsellis/data/IES_data/msas/alignments/aln/';
my $pathOut = '/home/dsellis/data/IES_data/msas/alignments/charMat/';
opendir(DH, $path) or die $!;
my @charMF = grep {/\.F\.dat$/} readdir(DH);

#read 2merge files and character matrices
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
	print $path.$mergeF,"\n";

# #print Matrix
	print "M:\n";
	print join("\t",@colNames),"\n";
	for(my $i = 0; $i <= $#M; $i++){ # new matrix has the same number of rows
	    print $rowNames[$i],"\t";
	    for(my $j = 0; $j <= $#{$M[0]}; $j++){
		print $M[$i][$j],"\t";
	    }
	    print "\n";
	}


	open IN, $path.$mergeF or die $!;
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
		    my $mergedStringS = join (',',@mergedString);
                    #new matrix gets merged column entries in place of the first of the merged ones
		    $NM[$i][$merged[0]] = $mergedStringS;
		}
	    }else{
		for(my $i = 0; $i <= $#M; $i++){
		    $NM[$i][$cols] = $M[$i][$cols];#new matrix has the same column
		}
	    }
	}
	close IN;
	open OUT, '>'.$pathOut.'cluster.'.$cluster.'.dat' or die $!;
	print "Cluster $cluster:\n";
	print OUT "geneName\t",join("\t",@NMcolNames),"\n";
	for(my $i = 0; $i <= $#M; $i++){ # new matrix has the same number of rows
	    print OUT $rowNames[$i],"\t";
	    for(my $j = 0; $j <= $#{$NM[0]}; $j++){
		if(!defined($NM[$i][$j])){
		    print $cluster,"\n";
		    use Data::Dumper;
		    print Dumper @NM;
		    die;
}
		print OUT $NM[$i][$j],"\t";
	    }
	    print OUT "\n";
	}
	close OUT;
    }
}
