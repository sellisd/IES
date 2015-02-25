#!/usr/bin/perl
use warnings;
use strict;

#correct start/end of floating IES
my $path = '/home/dsellis/data/IES_data/';

my @iesFiles = qw(
pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3
ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3
psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3
);
my $PtetMacF = $path.'ptetraurelia/Pte.ies.mac_seq';
my %PtmacSeq;
open PT, $PtetMacF or die $!;
while (my $line =<PT>){
    chomp $line;
    (my $iesID, my $mac_seq) = split "\t", $line;
    $PtmacSeq{$iesID} = $mac_seq;
}
close PT;
foreach my $file (@iesFiles){
    my $outputF = $file;
    $outputF =~ s/\.gff3/.fl.gff3/ or die;
    open OUT, '>'.$path.$outputF or die $!;
    open IN, $path.$file or die $!;
    while(my $line = <IN>){
	chomp $line;
	my @ar = split "\t", $line;
	my $start = $ar[3];
	my $end = $ar[4];
	my $annotation = $ar[8];
	my $seq;
	my $id;
	my $alt_seq;
	my $mac_seq;
	my @annotations = split ';', $annotation;
	my $counter = 0;
	foreach my $annot (@annotations){
	    if ($annot =~ /ID=(.*)/){
		$id = $1;
	    }elsif ($annot =~ /sequence=(.*)/){
		$seq = $1;
	    }elsif($annot =~ /alternative_IES_seq=(.*)/){
		$alt_seq = $1;
	    }elsif($annot =~ /mac_seq=(.*)/){
		$mac_seq = $1;	    
	    }
	}
	if(!defined($mac_seq)){
	    $mac_seq = $PtmacSeq{$id};
	}
	$mac_seq=~/([actg]+)([ACTG]+)[actg]+/ or die $id;
	my $window = length($1);
	my $ies_junction = length($2);
	my $macBorders = substr($mac_seq,$window-2, $ies_junction + 4);
	my $TAta; # bounded by TATA:
	my $taTA;
	if(substr($mac_seq, $window-2,2) eq 'ta'){
	    $taTA = 1;
	}
	if(substr($mac_seq, $window+$ies_junction,2) eq 'ta'){
	    $TAta = 1;
	}
	if($taTA){ #floating candidate
	    if(substr($seq,-4,4) eq 'TATA'){
		$start -= 2;
	    }
	}
	if($TAta){
	    if(substr($seq,0,4) eq 'TATA'){
		$end +=2;
	    }
	}
	$ar[3] = $start;
	$ar[4] = $end;
	my $newLine = join "\t", @ar;
	print OUT $newLine;
    }
    close IN;
}
	
