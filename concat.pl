#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
#use Bio::Align::Utilities qw(:all);
use lib'.';
use functions;

my $outFile = '/home/dsellis/data/IES/analysis/sgf/concat.fa';

my @species = qw/
Paramecium_primaurelia
Paramecium_biaurelia
Paramecium_tetraurelia
Paramecium_pentaurelia
Paramecium_sexaurelia
Paramecium_octaurelia
Paramecium_tredecaurelia
Paramecium_sonneborni
Paramecium_caudatum
Tetrahymena_thermophila
/;


my @end;
my @geneFamilies;
my $lengthSum = 0;
my %concat;
foreach my $inputFile (@ARGV){
    my $in = Bio::SeqIO->new(-file   => $inputFile,
			     -format => 'Fasta');
    my $length;
    my %aln;
    while(my $seqO = $in->next_seq()){
	my $species = gene2species($seqO->id);
	$aln{$species} = $seqO->seq();
	$length = $seqO->length();
    }
    foreach my $sp (@species){
	if(defined($aln{$sp})){
	    $concat{$sp} .= $aln{$sp};
	}else{
	    $concat{$sp} .= '-' x $length
	}
    }
    $lengthSum += $length;
    push @end, $lengthSum;
    push @geneFamilies, $inputFile;
}


# from end points recover start end coordinats bed format
my @start = @end;
pop @start; # remove last element
map($_++,@start);
unshift @start, 0; # adds the first start as 0

for(my $i = 0; $i <= $#start; $i++){
    printab($geneFamilies[$i], $start[$i], $end[$i]);
}

my $out = Bio::SeqIO->new(-file => '>'.$outFile,
			  -format => 'Fasta');
foreach my $sp (@species){
#    print $sp,' ', $concat{$sp},"\n";
    my $seqO = Bio::Seq->new(-id  => $sp,
 			     -seq => $concat{$sp});
    $out->write_seq($seqO);
}
