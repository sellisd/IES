#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use File::Path qw(make_path);
use Getopt::Long;
my $help;
my $home = '/home/dsellis/';
#default values
my $fastaOutPath = $home.'data/IES/analysis/msas/fasta/';
my $silixGroups = $home.'data/IES/analysis/allvsall/blastout/silix.output';
my $fastaAll = $home.'data/IES/analysis/protdb/allprot.fa';

my $usage = <<HERE;
Splits a fasta file to groups based on clustering from silix
Usage
./splitFasta.pl [OPTIONS]
where OPTIONS can be:
  out    : directory where output fasta files will be created
  silix  : path and filename of silix output file
  in     : path and filename of fasta files containing all sequences
  help   : this help screen

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'out=s' => \$fastaOutPath,
			      'silix=s' => \$silixGroups,
			      'in=s' =>\$fastaAll));
die $usage if $help;


#creates one fasta file for each silix group to be used for the multiple sequence alignments
#load groupings in memory

make_path($fastaOutPath) unless -d $fastaOutPath;

print "loading silix output in memory\n";
my %hash;

open IN, $silixGroups or die $!;
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
my $fastaFile = Bio::SeqIO->new(-file => $fastaAll,
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
