#!/usr/bin/perl
use warnings;
use strict;
my $cluster = 10189;
my $path = '/home/dsellis/data/IES_data/msas/alignments/aln/cluster.';
my $charMatF = $path.$cluster.'.F.dat';
#make character matrix based on nucleotide sequence coordinates
#    read character matrices (aa location and frame). if IES too close (calculate in bp) read IES sequence file and check ends of sequence and MAC sequence

open IN, $charMatF or die $!; #read character matrix
my $lineCounter = 0;
my @header;
my %geneIESH;
my %allNlocH; # hash to fit all unique locations in nucl. coordinates
while (my $line = <IN>){
    chomp $line;
    my @lineAr = split " ", $line;
    if($lineCounter == 0){ #header
	shift @lineAr; #remove genename
	@header = @lineAr;
# find if ies are closer than 2aa to each other
    }else{
# find if ies have a different insertion frame
	my $geneName = shift @lineAr;
	my @frames = @lineAr;
	my $colCounter = 0;
	foreach my $frame (@frames){
	    my $aaloc;
	    my $Nloc = 0;
	    if ($frame !=0){
		$aaloc = $header[$colCounter];
		$Nloc = $aaloc*3 - (3 - $frame);  #Nloc = aaloc*3 - (3-frame)	    
		$allNlocH{$Nloc} = 1;
                #create hash {geneName} = [nloc1, nloc2]
		if(defined()){
		    push @{$geneIESH{$geneName}}, $Nloc;
		}else{
		    $geneIESH{$geneName} = [$Nloc];
		}
	    }
	    $colCounter++;
	}
# find distances (in pb) of consecutive IES
# if they are too close to each other test for floating
    }
    $lineCounter++;
}
close IN;

# second pass go through again and print
my @newHeader;
#push @newHeader, 'geneName';
foreach my $nl (sort{$a <=> $b} keys %allNlocH){
    push @newHeader, $nl;
}

print 'geneName ',"@newHeader\n";
foreach my $gn (keys %geneIESH){
    print $gn,' ';
    foreach my $nloc (@newHeader){
	foreach my $iesL (@{$geneIESH{$gn}}){
	    if($iesL == $nloc){
#		die $iesL;
		print 1;
	    }else{
		print 0;
	    }
	    print ' ';
	}
    }
    print "\n";
}
 
