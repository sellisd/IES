#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Long;
use Bio::LocatableSeq;
my $help;
my @cds;
my $noterm;
my $usage = <<HERE;

read a protein multiple sequence alignment and replace the aligned proteins by nucleotide sequences
usage prot2nucl.pl [OPTIONS] prot1.aln.fa prot2.aln.fa
where OPTIONS can be:
    -cds:    fasta file with CDS. Multiple files can be specified with multiple -cds options
    -noterm: do not include termination codon, replace with gaps and remove all gap columns
    -help|?: this help screen

HERE

die $usage unless (GetOptions('help|?' => \$help,
			      'noterm' => \$noterm,
			      'cds=s'  => \@cds));
die $usage if $help;
die $usage if $#ARGV < 0;
my $prefixesR = initF();
my %seqIndexH;
print " build a hash of gene names and sequences\n";
foreach my $fileIN (@cds){
    print $fileIN, "\n";
    my $stream = Bio::SeqIO->new('-file'   => $fileIN,
				 '-format' => 'Fasta');
    while(my $seqO = $stream->next_seq){
	my $header = $seqO->display_id();
	my $sequence = $seqO->seq();
	#replace T in name with P
	$header = X2Y($header, $prefixesR, 'T', 'P');
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
    $outputF =~ s/\.aln\.fa$/.nucl.fa/ or die $!;
    my $alnStream = Bio::AlignIO->new('-file'   => $alnF,
				      '-format' => 'Fasta');
    my $newalnF = Bio::AlignIO->new('-file'   => '>'.$outputF,
				    '-format' => 'Fasta');
    my $newAlnO = Bio::SimpleAlign->new();
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
		    my $aa = substr($seq,$c,1);
		    if ($aa eq '-' or ($noterm and ($aa eq '*'))){
			$newseq .= '---';
		    }else{
			$newseq .= substr($nuclSeq,$curPos,3);
			$curPos+=3;
		    }
		}
		my $newLSO = Bio::LocatableSeq->new('-seq' => $newseq,
						    '-id' => $id);
		my $newSeqO = Bio::Seq->new('-display_id' => $id,
					    '-seq'        => $newseq);
		$newAlnO->add_seq($newLSO);
	    }else{
		die $id;
	    }
	}
	if($noterm){
	    my $noAllGapAlnO = $newAlnO->remove_columns(['all_gaps_columns']);
	    $newalnF->write_aln($noAllGapAlnO);
	}else{
	    $newalnF->write_aln($newAlnO);
	}
    }
}

