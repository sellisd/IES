#!/usr/bin/perl
use warnings;
use strict;

#prepare a catalog of IES names and DNA sequence of their junctions
# pte has a different gff3 file without mac-sequences, get those from Pte.ies.mac_seq

#read gff files
my @files = qw#
/home/dsellis/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3
/home/dsellis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3
/home/dsellis/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3 
#;

my %iesH;
my %ptetH;
open IN, '/home/dsellis/data/IES_data/ptetraurelia/Pte.ies.mac_seq' or die $!;
while (my $line = <IN>){
    chomp $line;
    (my $id, my $macseq) = split " ", $line;
    if(defined($ptetH{$id})){
	die "duplicate ids: $ptetH{$id}";
    }else{
	$ptetH{$id} = $macseq;
    }
}
close IN;

foreach my $file (@files){
    open IN, $file or die $!;
    while(my $line = <IN>){
	chomp $line;
	my $id;
	my $mac_seq;
	my $seq;
	my $str = (split " ", $line)[8];
	my @ar = split ';', $str;
	foreach my $pair(@ar){
	 if($pair =~ /ID=(.*)/){
	     $id = $1;
	 }elsif($pair =~ /mac_seq=(.*)/){
	     $mac_seq = $1;
	 }elsif($pair =~ /sequence=(.*)/){
	     $seq=$1;
	 }
	}
	if(defined($iesH{$id})){
	    die "duplicate IES ID: $iesH{$id}";
	}else{
	    if(!defined($mac_seq)){
		$mac_seq = $ptetH{$id}
	    }
	    $iesH{$id}=$mac_seq."\t".$seq;
	}
    }
    close IN;
}

foreach my $id (sort keys %iesH){
    print $id,"\t",$iesH{$id},"\n";
}
