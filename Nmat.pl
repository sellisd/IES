#!/usr/bin/perl
use warnings;
use strict;

# read character matrices in aminoacid coordinates and
# make character matrix in nucleotide coordinates with IES
# name in case of presence

my $path  = '/home/dsellis/data/IES_data/msas/alignments/aln/';

opendir(DH, $path) or die $!;
my @files = grep {/\.F\.dat/} readdir(DH);

foreach my $file(@files){
    $file =~ /cluster\.(\d+)\.F\.dat/;
    my $group = $1;
    my $fileWithIDs = 'cluster.'.$group.'.L.dat';
    my $outfile = 'cluster.'.$group.'.N.dat';
#    $outfile =~ s/\.F\.dat$/.N.dat/;
    print $outfile,"\n";
    my $charMatF = $path.$file;
#make character matrix based on nucleotide sequence coordinates
#    read character matrices (aa location and frame). if IES too close (calculate in bp) read IES sequence file and check ends of sequence and MAC sequence
    open IN, $charMatF or die $!; #read character matrix
    open IDF, $path.$fileWithIDs or die $!;
    my $lineCounter = 0;
    my @header;
    my %geneIESH;
    my %allNlocH; # hash to fit all unique locations in nucl. coordinates
    while (my $line = <IN>){
	my $lineID = <IDF>;
	chomp $line;
	chomp $lineID;
	my @lineAr = split " ", $line;
	my @lineArIDs = split " ", $lineID;
	if($lineCounter == 0){ #header
	    shift @lineAr; #remove genename
	    @header = @lineAr;
	}else{
	    my $geneName = shift @lineAr;
	    shift @lineArIDs; #keep arrays parallel
	    my @frames = @lineAr;
	    my $colCounter = 0;
	    foreach my $frame (@frames){
		my $aaloc;
		my $Nloc = 0;
		if ($frame != 0){
		    $aaloc = $header[$colCounter];
		    $Nloc = $aaloc*3 - (3 - $frame);  #Nloc = aaloc*3 - (3-frame) 
		    $allNlocH{$Nloc} = 1;
		    #create hash {geneName} = [nloc1, nloc2]
		    my $id = $lineArIDs[$colCounter];
		    if(defined($geneIESH{$geneName})){
			push @{$geneIESH{$geneName}}, {'ID' => $id, 'location'=>$Nloc};
		    }else{
			$geneIESH{$geneName} = [{'ID' => $id, 'location'=>$Nloc}];
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
    open OUT, '>'.$path.$outfile or die $!;
    print OUT 'geneName',"\t",join("\t",@newHeader),"\n";
    foreach my $gn (keys %geneIESH){
	print OUT $gn,"\t";
	foreach my $nloc (@newHeader){
	    my $flag = 0;
	    foreach my $iesL (@{$geneIESH{$gn}}){
		if($iesL->{'location'} == $nloc){
		    $flag = $iesL->{'ID'};
		}
	    }
	    print OUT $flag;
	    print OUT "\t";

	}
	print OUT "\n";
    }
    close OUT;
}
