#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

my $silixFile = '/Users/diamantis/data/IES_data/working/102N.silix';
my @proteinFiles = qw{
/Users/diamantis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.protein.fa
/Users/diamantis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.protein.fa
/Users/diamantis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.protein.fa
} ;

my @geneFiles = qw{
/Users/diamantis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gene.fa
/Users/diamantis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gene.fa
/Users/diamantis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.gene.fa
} ;

my $outputPFile = '>/Users/diamantis/data/IES_data/working/genes102P.fa';
my $outputNFile = '>/Users/diamantis/data/IES_data/working/genes102N.fa';
my $groupTarget = '15';
my %genes;
my %proteins; #same data but (?)faster access
#select a group from a silix clustering file and extract its protein sequence in a fasta file
open IN, $silixFile or die $!;

while(my $line = <IN>){
    chomp $line;
    (my $group, my $name) = split " ", $line;
    if ($group eq $groupTarget){
	$name =~ /P.*(P.*)T(\d+)/; #this depends on the output of extractExon102.pl if I correct it (for shorter output and less redundancy in the names I should adjust the regexp
	$genes{$1.'G'.$2} = 1;
	$proteins{$1.'P'.$2} = 1;
    }
}
close IN;

#use Data::Dumper;print Dumper %genes;die;
my $OUTP = Bio::SeqIO->new('-file' => $outputPFile,
			  '-format' => 'Fasta');
#read protein file(s)
foreach my $file (@proteinFiles){
    my $PF = Bio::SeqIO->new('-file' => $file,
			     '-format' => 'Fasta');
    while(my $seq = $PF->next_seq()){
	my $id = $seq->display_id();
#	print $id,' ', $seq->display_id();die;
	if (defined($proteins{$id})){
	    $OUTP->write_seq($seq);
	}
    }
}

#read gene file(s)

#use Data::Dumper;print Dumper %genes;die;
my $OUTN = Bio::SeqIO->new('-file' => $outputNFile,
			  '-format' => 'Fasta');
#read gene file(s)
foreach my $file (@geneFiles){
    my $NF = Bio::SeqIO->new('-file' => $file,
			     '-format' => 'Fasta');
    while(my $seq = $NF->next_seq()){
	my $id = $seq->display_id();
#	print $id,' ', $seq->display_id();die;
	if (defined($genes{$id})){
	    $OUTN->write_seq($seq);
	}
    }
}
