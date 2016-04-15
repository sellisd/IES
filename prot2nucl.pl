#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;

#read a protein multiple sequence alignment and replace the aligned proteins by nucleotide sequences

my @dnaF = qw#
/home/dsellis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.cds.fa
/home/dsellis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.cds.fa
/home/dsellis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.cds.fa
/home/dsellis/data/IES_data/pcaudatum_43c3d_annotation_v2.0/pcaudatum_43c3d_annotation_v2.0.cds.fa
/home/dsellis/data/IES_data/tthermophila/T_thermophila_June2014_CDS.fasta
#;

my %seqIndexH;
print " build a hash of gene names and sequences\n";
foreach my $fileIN (@dnaF){
    my $stream = Bio::SeqIO->new('-file'   => $fileIN,
				 '-format' => 'Fasta');
    while(my $seqO = $stream->next_seq){
	my $header = $seqO->display_id();
	my $sequence = $seqO->seq();
	#replace T in name with P
	$header =~ s/^(.+)T(\d+)$/$1P$2/;
	if(defined($seqIndexH{$header})){
	    die "duplicate gene name:", $seqIndexH{$header};
	    #make sure names are unique
	}else{
	    $seqIndexH{$header} = $sequence;
	}
    }
}

my @files = @ARGV;
print "read alignments\n";
foreach my $alnF (@files){
    next unless -f $alnF;
    die unless $alnF =~ /\.aln.fa$/;
    my $outputF = $alnF;
    print "  $outputF\n";
    $outputF =~ s/\.aln\.fa$/.nucl.fa/ or die $!;
    my $alnStream = Bio::AlignIO->new('-file'   => $alnF,
				      '-format' => 'Fasta');
    my $outputStream = Bio::SeqIO->new('-file'   => '>'.$outputF,
				       '-format' => 'Fasta');
    while(my $alnO = $alnStream->next_aln()){
	foreach my $seq ($alnO->each_seq()){
	    my $id = $seq->id();
	    my $seq = $seq->seq();
	    my $newseq = ''; # build here the aligned nucleotide sequence
	    my $curPos = 0; # pointer to the current location in nucleotide coordinates
	    if(defined($seqIndexH{$id})){	    
		my $nuclSeq = $seqIndexH{$id};
		# there should be no lowercase introns
		die if $nuclSeq =~ /[a-z]/; #drop any lowercase introns
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
					    '-seq'        => $newseq);
		$outputStream->write_seq($newSeqO);
	    }else{
		die;
	    }
	}
    }
}

