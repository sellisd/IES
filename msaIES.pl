#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
my %hash;
my @files;
my $home = '/home/dsellis/';
# Read alignment file and overlap file and print character matrices with IES locations and frames

my $dir = $ARGV[0];

opendir(DH, $dir) or die $!;
@files = grep { /\.nucl\.fa$/ } readdir(DH);

my @iesF = ($home.'data/IES_data/pbiaurelia/Pbi.IESinCDS',
	    $home.'data/IES_data/ptetraurelia/Pte.IESinCDS',
	    $home.'data/IES_data/psexaurelia/Pse.IESinCDS',
	    $home.'data/IES_data/pcaudatum_43c3d_annotation_v2.0/Pca.IESinCDS'
    );
foreach my $file (@iesF){  # load IES information in memory
    open IN, $file or die $!;
    while (my $line = <IN>){
    chomp $line;
    (my $gene, my $ies, my $start, my $end, my $length) = split " ",$line;
    $gene =~ /(.*)G(\d+)/;
    my $speciesAbr = $1;
    my $number = $2;
    $gene = $speciesAbr.'P'.$number; # replace with P for proteins 
    my $entry = {'ies' => $ies,
		 'species' => $speciesAbr,
		 'start'   => $start,
		 'end'   => $end,
		 'length'  => $length
    };
    if(defined($hash{$gene})){
	push @{$hash{$gene}},$entry;
    }else{
	$hash{$gene} =[$entry];
    }
    }
    close IN;
}


foreach my $alnF (sort @files){  # find IES coordinates in alignments
    my $alnIO = Bio::AlignIO->new(-file => $dir.$alnF,
				  -format=>'fasta'); # for the nucleotide alignments
    my %frameH;
    my @characterL;
    my @characterN;
    while(my $alnO = $alnIO->next_aln){
	my $seqNo = $alnO->num_sequences;
	print $alnF,"\t",$seqNo,"\t";
	print $alnO->percentage_identity(),"\t";
	#for each sequence in alignment
	my %charM;
	my %sortH; # for easy sorting by start
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $id = $seq->id();
	    if (defined($hash{$id})){      #find if it has IES
		foreach my $ies (@{$hash{$id}}){
		    my $start = $ies->{'start'};
		    my $end   = $ies->{'end'};
		    my $msaLocStart = $alnO->column_from_residue_number($id,$start);     #find the location of IES in the alignment
		    my $msaLocEnd   = $alnO->column_from_residue_number($id,$end);       #find the location of IES in the alignment
		    $charM{$msaLocStart.'.'.$msaLocEnd}{$id} = [$ies->{'length'}, $ies->{'ies'}];
		    $sortH{$msaLocStart.'.'.$msaLocEnd} = $msaLocStart; #sort by start
		}
	    }
	}

	# at least 1 IES present?
	my $colNo = scalar(keys %charM);
	if ($colNo>0){
	    print 1,"\n";
	}else{
	    print 0,"\n";
	}
	# second pass to print
	my $charMatrixFrameF = $alnF;
	my $charMatrixLengthF = $alnF;
	my $bedF = $alnF;
	$charMatrixFrameF =~ s/\.nucl\.fa/\.F.dat/ or die $!;
	$charMatrixLengthF =~ s/\.nucl\.fa/\.L.dat/ or die $!;
	$bedF =~ s/\.nucl\.fa/.bed/ or die $!;
	open OUTF, '>'.$dir.$charMatrixFrameF or die $!;
	open OUTL, '>'.$dir.$charMatrixLengthF or die $!;
	open BED, '>'.$dir.$bedF or die $!;
	print OUTF "geneName\t";
	print OUTL "geneName\t";
	$bedF =~ /cluster\.(\d+)\./;
	my $cluster = $1;
	my $colId = 0;
	foreach my $character (sort{$sortH{$a} <=> $sortH{$b}} keys %sortH){ # sort by start the IES locations
	    print OUTF $character,"\t";
	    print OUTL $character,"\t";
	    $character =~ /(\d+)\.(\d+)/;
	    my $start = $1;
	    my $end = $2;
	    print BED $cluster,"\t",$start,"\t", $end,"\t",$colId,"\n";
	    $colId++;
	}
	print OUTF "\n";
	print OUTL "\n";
	foreach my $seq ($alnO->each_seq()){
	    my $geneName = $seq->id();
	    print OUTF $geneName,"\t";
	    print OUTL $geneName,"\t";
	    my $id = $seq->id();
	    foreach my $state (sort {$sortH{$a} <=> $sortH{$b}} keys %sortH){
		if (defined($charM{$state}{$id})){
		    print OUTF ${$charM{$state}{$id}}[1],"\t"; # name
		    print OUTL ${$charM{$state}{$id}}[0],"\t"; # length
		}else{
		    if(substr($geneName,0,4) eq 'PBIA' or
		       substr($geneName,0,4) eq 'PTET' or
		       substr($geneName,0,4) eq 'PSEX' or
		       substr($geneName,0,4) eq 'PCAU'){
			print OUTF "0\t";
			print OUTL "0\t";
		    }else{ # or unknown
			print OUTF "NA\t";
			print OUTL "NA\t";		
		    }
		}
	    }
	    print OUTF "\n";
	    print OUTL "\n";
	}
	close OUTF;
	close OUTL;
	close BED;
    }
}
