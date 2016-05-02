#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use File::Spec::Functions qw(catfile);
use Bio::Tools::GFF;

my $tthP = '/home/dsellis/data/IES/thermophila/gene/';

#
# rename T. thermophila proteins
my $tthProt = 'T_thermophila_June2014.protein.fa';
my $tthCDS  = 'T_thermophila_June2014_CDS.fasta';
my $tthGene = 'T_thermophila_June2014.gene.fa';
my $tthGff  = 'T_thermophila_June2014.gff3';

my $outProt = 'tth.protein.fa';
my $outCds  = 'tth.cds.fa';
my $outGene = 'tth.gene.fa';
my $outGff  = 'tth.gff';

my %genes;
my $gffO = Bio::Tools::GFF->new('-file'        => catfile($tthP, $tthGff),
				'-gff_version' => 3);
while(my $feature = $gffO->next_feature()){
    my $scaffold = $feature->seq_id;
    my $type = $feature->primary_tag;
    if($type eq 'gene'){
	my $name = ($feature->get_tag_values('Name'))[0];
	my $note = ($feature->get_tag_values('Note'))[0];
	$genes{$name} = 1;
    }
}

# rename protein sequences
my $protIn = Bio::SeqIO->new('-file'   => catfile($tthP, $tthProt),
			     '-format' => 'Fasta');
my $protOut = Bio::SeqIO->new('-file'   => '>'.catfile($tthP, $outProt),
			      '-format' => 'Fasta');

while(my $protO = $protIn->next_seq){
    my $id =  $protO->primary_id;
    die if (!defined($genes{$id}));
    my $newSeqO = Bio::Seq->new(-seq => $protO->seq(),
			      -id  => $id,
	);
    $protOut->write_seq($newSeqO);
}

# rename CDS
my $cdsIn = Bio::SeqIO->new('-file'   => catfile($tthP, $tthCDS),
			    '-format' => 'Fasta');
my $cdsOut = Bio::SeqIO->new('-file'   => '>'.catfile($tthP, $outCds),
			     '-format' => 'Fasta');

while(my $cdsO = $cdsIn->next_seq){
    my $id =  $cdsO->primary_id;
    die if (!defined($genes{$id}));
    my $newSeqO = Bio::Seq->new(-seq => $cdsO->seq(),
				-id  => $id,
	);
    $cdsOut->write_seq($newSeqO);
}

# rename gene
my $geneIn = Bio::SeqIO->new('-file'   => catfile($tthP, $tthGene),
			     '-format' => 'Fasta');
my $geneOut = Bio::SeqIO->new('-file'   => '>'.catfile($tthP, $outGene),
			      '-format' => 'Fasta');

while(my $geneO = $geneIn->next_seq){
    my $id =  $geneO->primary_id;
    die if (!defined($genes{$id}));
    my $newSeqO = Bio::Seq->new(-seq => $geneO->seq(),
				-id  => $id,
	);
    $geneOut->write_seq($newSeqO);
}
