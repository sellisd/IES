#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

my $silixFile = '/Users/diamantis/data/IES_data/working/102P.silix';
my @proteinFiles = qw{
/Users/diamantis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.protein.fa
/Users/diamantis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.protein.fa
/Users/diamantis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.protein.fa
} ;

my $outputFile = '>/Users/diamantis/data/IES_data/working/genes102P.fa';
my $groupTarget = '1';
my %genes;
#select a group from a silix clustering file and extract its protein sequence in a fasta file
open IN, $silixFile or die $!;

while(my $line = <IN>){
    chomp $line;
    (my $group, my $name) = split " ", $line;
    if ($group eq $groupTarget){
	$name =~ /P.*(P.*)T(\d+)/; #this depends on the output of extractExon102.pl if I correct it (for shorter output and less redundancy in the names I should adjust the regexp
	$genes{$1.'P'.$2} = 1;
    }
}
close IN;

#use Data::Dumper;print Dumper %genes;die;
my $OUT = Bio::SeqIO->new('-file' => $outputFile,
			  '-format' => 'Fasta');
#read protein file(s)
foreach my $file (@proteinFiles){
    my $PF = Bio::SeqIO->new('-file' => $file,
			     '-format' => 'Fasta');
    while(my $seq = $PF->next_seq()){
	my $id = $seq->display_id();
#	print $id,' ', $seq->display_id();die;
	if (defined($genes{$id})){
	    $OUT->write_seq($seq);
	}
    }
}

#loop through sequences and print matching ones
