#!/usr/bin/perl
use warnings;
use strict;

# prepare a catalog of IES names and DNA sequence of their junctions ONLY for ies in proteins
# find strand for each gene
# find in which gene IES are located
# foreach gene find strand
# if strand is - revcomp IES and mac_seq

# pte has a different gff3 file without mac-sequences, get those from Pte.ies.mac_seq

# read gff3 files

my @anot = qw#
/home/dsellis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gff3
/home/dsellis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.gff3
/home/dsellis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gff3
#;

# read IESiCDS files
my @IESinCDS = qw#
/home/dsellis/data/IES_data/pbiaurelia/Pbi.IESinCDS
/home/dsellis/data/IES_data/psexaurelia/Pse.IESinCDS
/home/dsellis/data/IES_data/ptetraurelia/Pte.IESinCDS
#;

# read ies gff files
my @files = qw#
/home/dsellis/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3
/home/dsellis/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3
/home/dsellis/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3 
#;

my %iesH;        # {iesid  => macseq\tseq
my %ptetH;       # {iesid  => macseq} only for Pte $#!@
my %iesincdsH;   # {iesid  => geneid}
my %geneStrandH; # {geneid => strand}
#read genes with IES
foreach my $file (@IESinCDS){
    open IN, $file or die $!;
    while(my $line = <IN>){
	chomp $line;
	(my $gene, my $ies) = (split " ", $line)[0,1];
	$iesincdsH{$ies} = $gene;
    }
    close IN;
}

#read gene strand
foreach my $file (@anot){
    open IN, $file or die $!;
    while(my $line = <IN>){
	next if substr($line,0,1) eq '#';
	chomp $line;
	(my $feature, my $strand, my $str) = (split " ", $line)[2,6,8];
	next unless $feature eq 'gene';
	my @ar = split ';', $str;
	my $id;
	foreach my $pair (@ar){
	    if(substr($pair,0,3) eq 'ID='){
		$id = substr($pair,3);	
	    }
	}
	$geneStrandH{$id} = $strand;
    }
    close IN;
}

#read Ptetraurelia IES sequences
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

#read IES
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
	    if(defined($iesincdsH{$id})){ #if IES is in a gene
		my $strand = $geneStrandH{$iesincdsH{$id}};
		if($strand eq '+'){
		    # do nothing
		}elsif($strand eq '-'){
		    #revcomp seq and mac_seq
		    print $id,' ',$mac_seq,' ',$seq,"\n";
		    $mac_seq =~ tr/ATGCatgc/TACGtacg/;
		    $seq =~  tr/ATGCatgc/TACGtacg/;
		    $mac_seq = reverse($mac_seq);
		    $seq = reverse($seq);
		}else{
		    die;
		}
		$iesH{$id} = $mac_seq."\t".$seq; 
	    }
	}
    }
    close IN;
}


foreach my $id (sort keys %iesH){
    print $id,"\t",$iesH{$id},"\n";
}
