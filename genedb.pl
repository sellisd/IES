#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Getopt::Long;
use Bio::Tools::GFF;

my $help;
my $gff;
my $silixOutF;
my $origff;
my $usage = <<HERE;

combine the result of analyses to create the gene table for iesDB
usage genedb.pl [OPTIONS]
where OPTIONS can be:
    -silixout : file with the silix output containing the gene family assignments
    -gff :      gff file with gene annotation
    -origff :   original (unfiltered) gff file with gene annotation
    -help|?: this help screen

HERE
    
die $usage unless (GetOptions('help|?'      => \$help,
			      'silixout=s'  => \$silixOutF,
			      'origff=s'    => \$origff,
 			      'gff=s'       => \$gff,
		   ));

die $usage if $help;
# read silix output and get genes per cluster
open S, $silixOutF or die $!;
my %gf; #gene families
my $prefixesR = initF();

while(my $line = <S>){
    (my $geneFamily, my $protName) = split " ", $line;
    my $geneName = prot2gene($protName, $prefixesR);
    $gf{$geneName} = $geneFamily;
}
close S;

# prepare gene info db
#with output information
#id(unique) protName geneFamily species scaffold begin end filtered

# find out filtered genes
my $gffIO = Bio::Tools::GFF->new('-file' => $gff,
				 '-gff_version' => 3);
my %passedFilter; #genes that passed the scaffold filtering step
while(my $feature = $gffIO->next_feature()){
    next unless $feature->primary_tag eq 'gene';
    my $id = $feature->primary_id;
    die if defined($passedFilter{$id});
    $passedFilter{$id} = 1;
}

# read gff file
my $gffI = Bio::Tools::GFF->new('-file' => $origff,
 				'-gff_version' => 3);
my %genesdb;
printab('id', 'protName', 'geneFamily', 'scaffold', 'start', 'end', 'passedFilter');
while(my $feature = $gffI->next_feature()){
    my $scaffold = $feature->seq_id;
    my $start = $feature->start;
    my $end = $feature->end;
    my $type = $feature->primary_tag;
    my $id = $feature->primary_id;
    next unless $type eq 'gene';
    # foreach my $tag ($feature->get_all_tags()){
    # 	my @values = $feature->get_tag_values($tag);
    # }
    printab($id, X2Y($id, $prefixesR, 'G', 'P'), (defined($gf{$id})?$gf{$id}:'NA'), $scaffold, $start, $end, (defined($passedFilter{$id})?1:0));
}
