#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
# find length of sequences and print to database
# calculate coverage from blast file

my $iesF = '/home/dsellis/data/IES_data/ptetraurelia/ies.fa';
my $seqIO = Bio::SeqIO->new('-file' => $iesF,
			   '-format' => 'fasta');
my %lengthsH;
while(my $seqO = $seqIO->next_seq){
    my $id = $seqO->display_id;
    my $length = $seqO->length;
    $lengthsH{$id} = $length;
}
open LDB, '>/home/dsellis/data/IES_data/iesnet/iesLengths.tab' or die $!;
foreach my $ies (sort keys %lengthsH){
    print LDB $ies,"\t",$lengthsH{$ies},"\n";
}
close LDB;

my $filteredF = '/home/dsellis/data/IES_data/iesnet/ptet.filtered';
my %covM; # hash of coverage matrices {ies1 => [0,1,1,1],
          #                            ies2 => [1,1,1,2]}
open IN, $filteredF or die $!;
while(my $line = <IN>){ # read blast output line - by line
    chomp $line;
    (my $qseqid, my $sseqid, my $pident, my $length, my $mismatch, my $gapopen, my $qstart, my $qend, my $sstart, my $send, my $evalue, my $bitscore) = split " ", $line; # parse columns
    if(!defined($covM{$qseqid})){# the firt time we find an IES create an empty matrix
	my @arr = (0) x $lengthsH{$qseqid};
	$covM{$qseqid} = \@arr;
    }
    #add the coverage
    for(my $i = $qstart-1; $i<=$qend-1; $i++){ #blast nucleotide coordinates are 1-based while the matrix is 0-based
	@{$covM{$qseqid}}[$i]++;
    }

}
close IN;
 # use Data::Dumper;
 # print "@{$covM{'IESPGM.PTET51.1.1.10090'}}\n";
 # die;
my $covX = 0;
while(1){
    if($covX == 0){ #print header with IES names
	print join("\t",sort(keys(%covM)));
#	foreach my $cm (keys %covM){
#	    print $cm,"\t";
#	}
	print "\n";
    }else{
# print covX percent
	my @row;
	my $flag = 0;
	foreach my $cm (sort keys %covM){
	    my $length = $#{$covM{$cm}} + 1; # length of nucleotides = index of last element + 1
	    my $sum = 0;
	    foreach my $element (@{$covM{$cm}}){
		if($element >= $covX){
		    $sum++;	
		}
	    }
	    if ($sum > 0){
		$flag = 1;
	    }
	    my $perc = $sum/$length;
	    push @row,$perc;
	}
	if($flag == 1){
	    print join"\t",@row;
	    print "\n";
	}else{
	    print "\n";
	    last;
	}
#    print $#row,"\n";
    }
    $covX++;
 
}
