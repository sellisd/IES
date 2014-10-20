#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
#create multiple sequence alignments of silix groups

my $dataPath = '/Users/diamantis/data/IES_data/';

#load groupings in memory
my %hash;
my $silixOutput = $dataPath.'working/silix.output';
print $silixOutput;
open IN, $silixOutput or die $!;
while (my $line = <IN>){
    chomp $line;
    (my $group, my $id) = split " ", $line;
    if(defined($hash{$group})){
	push @{$hash{$group}}, $id;
    }else{
	$hash{$group}=[$id];
    }
}
close IN;

#load all protein sequences in memory
my %proteins;
my @speciesFiles = ('pbiaurelia/biaurelia_V1-4_annotation_v1.protein.fa','ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.4.protein.fa','psexaurelia/sexaurelia_AZ8-4_annotation_v1.protein.fa');
foreach my $inP (@speciesFiles){
    my $file = Bio::SeqIO->new('-file'=>$dataPath.$inP,
			       '-format' => 'Fasta');
    while(my $seqO = $file->next_seq()){
	my $header = $seqO->primary_id();
	my $sequence = $seqO->seq();
	if(defined($proteins{$header})){
	    die "synonymia?";
	}else{
	    $proteins{$header}  = $sequence;
	}
    }
}
use Data::Dumper;

#foreach group make a fasta file with the sequencies
#make list of groups with more than one sequence and align only those
foreach my $group (keys %hash){
    open OUT, '>'.$dataPath.'msas/group.'.$group.'.fa' or die $!;
    foreach my $proteinId (@{$hash{$group}}){
	if(defined($proteins{$proteinId})){
     	    print OUT ">$proteinId\n$proteins{$proteinId}\n";
     	}else{
    	    die "$group empty group?";
     	}
    }
    close OUT;
}

#in R
#a<-read.table("working/silix.output")
#read silix output line by line
#make hash with group and proteins
#read protein files and make one fasta file for each group
#run t-coffee for each fasta file
