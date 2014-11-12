#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use File::Path qw(make_path);

my $dataPath = '/Users/diamantis/data/IES_data/';
my $fastaOutPath = '/Users/diamantis/data/IES_data/msas/fasta/';

#creates one fasta file for each silix group to be used for the multiple sequence alignments
#load groupings in memory

make_path($fastaOutPath) unless -d $fastaOutPath;

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

# read combined.fa
my %seqH;
print "loading protein sequencies in memory\n";
my $fastaFile = Bio::SeqIO->new(-file => $dataPath.'working/combined.fa',
				-format => 'Fasta');
while (my $seqO = $fastaFile->next_seq){
    my $name = $seqO->display_id;
    if(defined($seqH{$name})){
	die "Error: names not unique: $name";
    }else{
	$seqH{$name} = $seqO;
    }
}

# for each family print proteins, but only if family has 2 or more members
foreach my $cluster (keys %hash){
    if($#{$hash{$cluster}}>1){
	my $fastaOutF = Bio::SeqIO->new(-file => '>'.$fastaOutPath.'cluster.'.$cluster.'.fa',
					-format => 'Fasta');
	foreach my $protein (@{$hash{$cluster}}){
	    $fastaOutF->write_seq($seqH{$protein});
	}
    }
}
