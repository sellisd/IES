#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;

#Alternative approach: alternatively from gff3 get all 102CDS from annotation get their sequence and translate
my $home = '/home/dsellis/';
my @species = qw/Pbi Pte Pse/;
my $PF = Bio::SeqIO->new('-file' => '>'.$home.'data/IES_data/working/exon102P.fa',
			   '-format' => 'fasta');
my $NF = Bio::SeqIO->new('-file' => '>'.$home.'data/IES_data/working/exon102N.fa',
			   '-format' => 'fasta');
my $tabOutF = '>'.$home.'data/IES_data/working/exon102.tab';

#find the 102 exons and print a summary file with
#id start stop length coding frame
#and two fasta files one with nucleotide and the aminoacid sequence
#open TO, $tabOutF or die $!;
my %scaffoldH;
foreach my $species (@species){
    my $dataPath = $home.'data/IES_data/';
    print $species,"\n";
    my $gff3;
    my $scaffoldsF;
    if($species eq 'Pbi'){
	$dataPath = $dataPath.'pbiaurelia/';
	$gff3 = $dataPath.'pbiaurelia_V1-4_annotation_v2.0.gff3';
	$scaffoldsF = $dataPath.'biaurelia_V1-4_assembly_v1.fasta';
    }elsif($species eq 'Pte'){
	$dataPath = $dataPath.'ptetraurelia/';
	$gff3 = $dataPath.'ptetraurelia_mac_51_annotation_v2.0.gff3';
	$scaffoldsF = $dataPath.'ptetraurelia_mac_51.fa';
    }elsif($species eq 'Pse'){
	$dataPath = $dataPath.'psexaurelia/';
	$gff3 = $dataPath.'psexaurelia_AZ8-4_annotation_v2.0.gff3';
	$scaffoldsF = $dataPath.'sexaurelia_AZ8-4_assembly_v1.fasta';
    }else{
	die;
    }
    
    #open sequence files and fill hashes
    print "read scaffolds from $scaffoldsF\n";
    my $scaffoldIn = Bio::SeqIO->new('-file' => $scaffoldsF,
				     '-format' => 'fasta');
    while(my $scaffoldSeq = $scaffoldIn->next_seq){
	my $scId = $scaffoldSeq->display_id;
	$scId =~ /(scaffold_?[_\d]+)/ or die $scId; #Everything is in MAC coordinates (without IES) scaffold names are not consistent across species
	$scaffoldH{$1} = $scaffoldSeq;
    }

    my $gff3In = Bio::Tools::GFF->new('-file' => $gff3,
				      '-gff_version' => 3);
    
    print "parse gff3 features\n";
#read gff3 file and build genbank entries
    while(my $feature = $gff3In->next_feature()){ # one line at a time
	if ($feature->primary_tag() eq 'CDS'){
	    if($feature->length() == 102){
		my $scaffold = $feature->seq_id();
		my $start = $feature->start();
		my $end = $feature->end();
		my $strand = $feature->strand();
		my $phase = $feature->frame();
		my $gene = ($feature->get_tag_values('Parent'))[0];
		my $id = $species.'_'.$scaffold.'_'.$gene.'_'.$start.'_'.$end;
		my $newNSeq = Bio::Seq->new('-seq' => $scaffoldH{$scaffold}->subseq($start,$end),
					    '-id' => $id,
					    '-alphabet' => 'dna'
		    );
		my $newPSeq;
		if($strand == 1){
		    $newPSeq = $newNSeq->translate(-frame=>$phase,
						   -codontable_id => 6);
		}elsif($strand == -1){
		    $newPSeq = $newNSeq->revcom()->translate(-frame=>$phase,
							     -codontable_id => 6);
		}else{
		    die;
		}
		$newPSeq->id($id);
		$newPSeq->alphabet('protein');	
		$PF->write_seq($newPSeq);
		$NF->write_seq($newNSeq);
	    }
	    
	}
    }

}
   
#     my $geneIn = Bio::SeqIO->new('-file' => $gene,
# 				 '-format' => 'genBank');
#     while(my $seqO = $geneIn->next_seq){
# 	my $scaffold = $seqO->accession_number();
# 	my $species = $seqO->species();
# 	foreach my $featureO ($seqO->get_SeqFeatures()){
# 	    if($featureO->primary_tag() eq 'CDS'){


# 		my $length = $featureO->length();
# 		my $start = $featureO->start();
# 		my $end = $featureO->end();
# 		my $gene = ($featureO->get_tag_values('gene'))[0];
# 		my $codon_start = ($featureO->get_tag_values('codon_start'))[0];
# 		my $id = 'P_'.$species->species().'_'.$seqO->accession_number().'_'.$gene.'_'.$start.'_'.$end;
# 		if($length == 102){
# 		    my $newNSeq = Bio::Seq->new('-seq' => $seqO->subseq($start,$end),
# 						'-id' => $id,
# 						'-alphabet' => 'dna'
# 			);
# 		    my $newPSeq = $newNSeq->translate(-frame=>$codon_start,
# 							 -codontable_id => 6);
# 		    $newPSeq->id($id);
# 		    $newPSeq->alphabet('protein');
# 		    $PF->write_seq($newPSeq);
# 		    $NF->write_seq($newNSeq);
# 		    print TO 'P.'.$species->species()."\t$scaffold\t$gene\t$start\t$end\t$length\t$codon_start\n";
# 		}
# 	    }
# 	}
#     }
# }
# close TO;

