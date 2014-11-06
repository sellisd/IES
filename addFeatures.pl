#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Species;
use Getopt::Long;
my $help;
my $species3abr;
my $usage = <<HERE;

add features to genbank files
usage addFeatures.pl -species Pab
HERE

# parse Pab.IES.gbkn
#       silix.output
#       Pab.overlap

# output Pib.seq CDS /gene_family = X
#                IES_junction = in_gene = GeneName

#parse Pab.overlap
die $usage unless (GetOptions('help|?' => \$help,
			      'species=s' => \$species3abr));
die $usage if $help;
my $dataPath = '/Users/diamantis/data/IES_data/';
my $silixOutF = '/Users/diamantis/data/IES_data/working/silix.output';

if ($species3abr eq 'Pbi'){
    $dataPath = $dataPath.'pbiaurelia/';
}elsif($species3abr eq 'Pte'){
    $dataPath = $dataPath.'ptetraurelia/'; 
}elsif($species3abr eq 'Pse'){
    $dataPath = $dataPath.'psexaurelia/';
}else{
    die;
}
my $overlapF = $species3abr.'.overlap';
my $genbankF = $species3abr.'.IES.gnbk';
my $gnbkOut = $species3abr.'.seq';
my $iesOutF = $species3abr.'.ies.seq';
my $iesF = $species3abr.'.ies';

my %iesH;
my %geneH;
print "parse ies overlap\n";
open OV, $dataPath.$overlapF or die $!;
while(my $line = <OV>){
    (my $gene, my $ies) = (split " ", $line)[3,7];
    $iesH{$ies} = $gene;
}
close OV;

#parse silix.output
print "parse silix output\n";
open SO, $silixOutF or die $!;
while(my $line = <SO>){
    (my $group, my $gene) =  split " ", $line;
#    from protein to gene
    $gene =~ s/(.*\.)P(\d+)/$1G$2/;
    $geneH{$gene} = $group;
}
close SO;

#parse genbank files
my $gnbkIn = Bio::SeqIO->new('-file' => $dataPath.$genbankF,
			     '-format' => 'genbank');
my $iesIn = Bio::SeqIO->new('-file' => $dataPath.$iesF,
			    '-format' => 'genbank');
my $nbkOut = Bio::SeqIO->new('-file' => '>'.$dataPath.$gnbkOut,
			     '-format' => 'genbank');
my $iesOut = Bio::SeqIO->new('-file' => '>'.$dataPath.$iesOutF,
			     '-format' => 'genbank');

print "parse genbank files\n";
while(my $seqO = $gnbkIn -> next_seq()){
    my $scaffold = $seqO->accession_number();
    my $species = $seqO->species();
    print "  scaffold ", $scaffold,"\n";
    foreach my $featureO ($seqO->get_SeqFeatures()){
	if($featureO->primary_tag() eq 'gene'){
	    my $geneName = ($featureO->get_tag_values('gene'))[0];
	    if(defined($geneH{$geneName})){
		$featureO->add_tag_value('gene_family' => $geneH{$geneName}); #add silix grouping
	    }
	}
	if($featureO->primary_tag() eq 'IES_junction'){ #add gene in which it is located
	    my $iesName = ($featureO->get_tag_values('id'))[0];
	    if(defined($iesH{$iesName})){
		$featureO->add_tag_value('in_gene' => $iesH{$iesName});
	    }
	}
    }
    $nbkOut->write_seq($seqO);
}

while(my $seqO = $iesIn -> next_seq()){
    foreach my $featureO ($seqO->get_SeqFeatures()){
	if($featureO->primary_tag() eq 'IES'){
	    my $iesName = ($featureO->get_tag_values('id'))[0];
	    if(defined($iesH{$iesName})){
		$featureO->add_tag_value('in_gene' => $iesH{$iesName});
	    }
	}
    }
    $iesOut->write_seq($seqO);
}

