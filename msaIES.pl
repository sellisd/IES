#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
my %hash;
my @files;
my $home = '/home/dsellis/';
#read alignment file and overlap file and print character matrices with IES locations and frames

my $dir = $ARGV[0];

opendir(DH, $dir) or die $!;
@files = grep { /nucl.fa/ } readdir(DH);

my @iesF = ($home.'data/IES_data/pbiaurelia/Pbi.IESinCDS',
	    $home.'data/IES_data/ptetraurelia/Pte.IESinCDS',
	    $home.'data/IES_data/psexaurelia/Pse.IESinCDS');
foreach my $file (@iesF){
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


foreach my $alnF (sort @files){
    my $alnIO = Bio::AlignIO->new(-file => $dir.$alnF,
				  -format=>'fasta'); # for the nucleotide alignments
#   while(my $alnO = $alnIO->next_aln){
# 	my $seqNo = $alnO->num_sequences;
# 	print $alnF,' ',$seqNo,"\t";
# 	print $alnO->percentage_identity(),"\n";
# 	foreach my $seq ($alnO->each_seq()) {     #find which genes
# 	    my $id = $seq->id();
# 	    print '   ', $id,"\n";
# 	}
#    }
# }
# die;
# if(0){my $alnIO;my $alnO;my $alnF;
    my %frameH;
    my @characterL;
    my @characterN;
    while(my $alnO = $alnIO->next_aln){
	my $seqNo = $alnO->num_sequences;
	print $alnF,' ',$seqNo,"\t";
	print $alnO->percentage_identity(),"\n";
	next if ($seqNo < 3);
	next if($alnO->percentage_identity() < 80);
	#for each sequence in alignment
	my %charM;
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $id = $seq->id();
#	    $id =~s/P.._scaffold.*?(P.*)T(\d+)_.*/$1G$2/;
	    if (defined($hash{$id})){      #find if it has IES
		foreach my $ies (@{$hash{$id}}){
#		    my $aaLoc = $ies->{'aaLoc'}; 
##		    my $ntLoc = $ies->{'start'};
		    my $start = $ies->{'start'};
		    my $end   = $ies->{'end'};
		    my $msaLocStart = $alnO->column_from_residue_number($id,$start);     #find the location of IES in the alignment
		    my $msaLocEnd   = $alnO->column_from_residue_number($id,$end);     #find the location of IES in the alignment
		    #	  print $msaLoc,'(',$ies->{'frame'},', ',$aaLoc,', ',$ies->{'ies'},")\t";
	#	    print $msaLoc,"\t";
#		    $charM{$msaLoc}{$id} = [$ies->{'frame'}, $ies->{'length'}, $ies->{'ies'}];
		    $charM{$msaLocStart.'.'.$msaLocEnd}{$id} = [$ies->{'length'}, $ies->{'ies'}];
		}
	    }
	 #   print "\n";
	}
	#only print if at least 1 IES is present
	next unless ((scalar(keys %charM)) > 0);
	#print "IES presence character matrix \n";
	#second pass to print
	my $charMatrixFrameF = $alnF;
	my $charMatrixLengthF = $alnF;
	$charMatrixFrameF =~ s/\.aln/\.F.dat/;
	$charMatrixLengthF =~ s/\.aln/\.L.dat/;
	open OUTF, '>'.$dir.$charMatrixFrameF or die $!;
	open OUTL, '>'.$dir.$charMatrixLengthF or die $!;
	print OUTF "geneName\t";
	print OUTL "geneName\t";
#	foreach my $geneName (sort {$a<=>$b} keys %charM){
	foreach my $character (sort keys %charM){ # unsorted
	    print OUTF $character,"\t";
	    print OUTL $character,"\t";
	}
	print OUTF "\n";
	print OUTL "\n";
	foreach my $seq ($alnO->each_seq()){
	    print OUTF $seq->id(),"\t";
	    print OUTL $seq->id(),"\t";
	    my $id = $seq->id();
#	    foreach my $state (sort {$a<=>$b} keys %charM){
	    foreach my $state (sort keys %charM){
#	    print $charM{$state},"\n";
		if (defined($charM{$state}{$id})){
		    print OUTF ${$charM{$state}{$id}}[1],"\t"; # length
		    print OUTL ${$charM{$state}{$id}}[1],"\t"; # name
		}else{
		    print OUTF "0\t";
		    print OUTL "0\t";
		}
	    }
	    print OUTF "\n";
	    print OUTL "\n";
	}
	close OUTF;
	close OUTL;
    }
}
