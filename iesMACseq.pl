#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
# find location of IES on assembled genome and extract the surrounding MAC region

#read fasta file
#read ies locations
my $file = 'internal_eliminated_sequence_PGM_IES51.pt_51.gff3';
my $path = '/home/dsellis/data/IES_data/ptetraurelia/';
my $outputF = 'Pte.ies.mac_seq';
my $window = 15;
my %iesH;
print "reading IES information\n";
open IN, $path.$file or die $!;
while (my $line = <IN>){
	chomp $line;
	my @ar = split "\t", $line;
	my $scaffold = $ar[0];
	my $annotation = $ar[8];
	my $id;
	my $start = $ar[3];
	my $end = $ar[4];
	my @annotations = split ';', $annotation;
	foreach my $annot (@annotations){
	    if ($annot =~ /ID=(.*)/){
		$id = $1;
	    }
	}
	$iesH{$scaffold}{$id} = {
	    'start'    => $start,
	    'end'      => $end
	};
}
close IN;

print "reading Fasta\n";
my $FA = Bio::SeqIO->new('-file' => $path.'ptetraurelia_mac_51.fa',
			 '-format' => 'fasta');
print "processing and printing\n";
open OUT, '>'.$path.$outputF or die $!;
while(my $seqO = $FA->next_seq){
    my $scID = $seqO->display_id;
    if(defined($iesH{$scID})){
	foreach my $iesID (keys %{$iesH{$scID}}){
	    my $start       = $iesH{$scID}{$iesID}->{'start'};
	    my $end         = $iesH{$scID}{$iesID}->{'end'};
	    my $iesJunction = uc($seqO->subseq($start,$end));
	    my $upstream    = lc($seqO->subseq($start - $window, $start-1));
	    my $downstream  = lc($seqO->subseq($end+1, $end + $window));
	    print OUT $iesID,"\t",$upstream.$iesJunction.$downstream,"\n";
	}

    }
}

close OUT;
