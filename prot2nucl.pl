#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
#read a protein multiple sequence alignment and replace the aligned proteins by nucleotide sequences

my @dnaF = qw#
/home/dsellis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gene.fa
/home/dsellis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gene.fa
/home/dsellis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.gene.fa
/home/dsellis/data/IES_data/pcaudatum_43c3d_annotation_v2.0/pcaudatum_43c3d_annotation_v2.0.gene.fa
#;

my $alnF = '/home/dsellis/data/IES_data/msas/alignments/aln/cluster.66.aln';
my $outputF = $alnF;
$outputF =~ s/\.aln/.nucl.fa/;
my %seqIndexH;
#build a hash of gene names and sequences
foreach my $fileIN (@dnaF){
    my $stream = Bio::SeqIO->new('-file'   => $fileIN,
				 '-format' => 'Fasta');
    while(my $seqO = $stream->next_seq){
	my $header = $seqO->display_id();
	my $sequence = $seqO->seq();
        #replace G in name with P
	$header =~ s/^(.+)G(\d+)$/$1P$2/;
	$seqIndexH{$header} = $sequence;
    }
}

#READ msa
#replace
#	use Data::Dumper; print Dumper %seqIndexH;
my $alnStream = Bio::AlignIO->new('-file'   => $alnF,
				  '-format' => 'clustalw');
my $outputStream = Bio::SeqIO->new('-file'   => '>'.$outputF,
				   '-format' => 'Fasta');
#my $outputStream = Bio::AlignIO->new('-file'   => '>'.$outputF,
#			       '-format' => 'clustalw');
while(my $alnO = $alnStream->next_aln()){
#    my @newSeqObjects;
    foreach my $seq ($alnO->each_seq()){
	my $id = $seq->id();
	my $seq = $seq->seq();
	my $newseq = ''; # build here the aligned nucleotide sequence
	my $curPos = 0; # pointer to the current location in nucleotide coordinates
	if(defined($seqIndexH{$id})){	    
	    my $nuclSeq = $seqIndexH{$id};
	    $nuclSeq =~ s/[a-z]+//g; #drop any lowercase introns
	    # for each aminoacid (unless it is a gap) print three nucleotides
	    for (my $c = 0; $c < length($seq); $c++){
		if (substr($seq,$c,1) eq '-'){
		    $newseq .= '---';
		}else{
		    $newseq .= substr($nuclSeq,$curPos,3);
		    $curPos+=3;
		}
	    }
	    my $newSeqO = Bio::Seq->new('-display_id' => $id,
					'-seq'              => $newseq);
#	    push @newSeqObjects, $newSeqO;
	    $outputStream->write_seq($newSeqO);
	}
    }
#    my $newAlnO = Bio::SimpleAlign->new( '-source' => 'prot2nucl.pl',
#					 '-seqs'   => \@newSeqObjects);
}
#$outputStream->write_seq($newAln);
