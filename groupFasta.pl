#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
my $home = '/home/dsellis/';
my $dataPath = $home.'data/IES_data/';
my $fastaOutPath = $home.'data/IES_data/msas/fasta/';

#make required directories
mkdir $fastaOutPath unless (-d $fastaOutPath);

#load groupings in memory
print "loading silix output in memory\n";
my %hash;
my $silixOutput = $dataPath.'working/silix.output';

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
print "loading protein sequences in memory\n";
my %proteins;
my @speciesFiles = ('pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.protein.fa','ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.protein.fa','psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.protein.fa');
foreach my $inP (@speciesFiles){
    print "   from $inP\n";
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

#foreach group make a fasta file with the sequencies
#make list of groups with more than one sequence and align only those
my $groupNo = keys %hash;
print "making fasta files\n";
my @toAlign;
foreach my $group (keys %hash){
    if($#{$hash{$group}}>1){ #make fasta files and align only if thee or more than one sequence in a group
	push @toAlign, $group;
	open OUT, '>'.$fastaOutPath.'/group.'.$group.'.fa' or die $!;
	foreach my $proteinId (@{$hash{$group}}){
	    if(defined($proteins{$proteinId})){
		print OUT ">$proteinId\n$proteins{$proteinId}\n";
	    }else{
		die "$group empty group?";
	    }
	}
	close OUT;
    }
}
