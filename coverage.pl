#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
# find length of sequences
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

my $filteredF = '/home/dsellis/data/IES_data/iesnet/ptet.filtered';
my %covM; # hash of coverage matrices
open IN, $filteredF or die $!;
while(my $line = <IN>){
    chomp $line;
    (my $qseqid, my $sseqid, my $pident, my $length, my $mismatch, my $gapopen, my $qstart, my $qend, my $sstart, my $send, my $evalue, my $bitscore) = split " ", $line;
    if(defined($covM{$qseqid})){
	for(my $i = $qstart; $i<=$qend; $i++){
	    @{$covM{$qseqid}}[$i]++;
	}
    }else{
	my @arr = (0) x $lengthsH{$qseqid};
	$covM{$qseqid} = \@arr;
    }
}
close IN;
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
	    my $length = $#{$covM{$cm}};
	    my $sum = 0;
	    foreach my $element (@{$covM{$cm}}){
		if($element == $covX){
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
