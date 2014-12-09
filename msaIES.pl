#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
my %hash;
my @files;
my $home = '/home/dsellis/';
#read alignment file and overlap file and print character matrices with IES locations and frames

my $dir = $ARGV[0]; #start with one alignment file, then modify to use all alignment files in a directory
opendir(DH, $dir) or die $!;
@files = grep { /aln/ } readdir(DH);
#print @files," \n";die;
my @iesF = ($home.'data/IES_data/pbiaurelia/Pbi.IESinCDS',
	    $home.'data/IES_data/ptetraurelia/Pte.IESinCDS',
	    $home.'data/IES_data/psexaurelia/Pse.IESinCDS');
foreach my $file (@iesF){
  open IN, $file or die $!;
  while (my $line = <IN>){
    chomp $line;
    (my $gene, my $ies, my $aaLoc, my $frame) = split " ",$line;
    $gene =~ /(.*)G(\d+)/;
    my $speciesAbr = $1;
    my $number = $2;
    $gene = $speciesAbr.'P'.$number; # replace with P for proteins 
    my $entry = {'ies' => $ies,
		 'species' => $speciesAbr,
		 'aaLoc' => $aaLoc,
		 'frame' =>$frame,
		};
      if(defined($hash{$gene})){
      push @{$hash{$gene}},$entry;
    }else{
      $hash{$gene} =[$entry];
    }
  }
  close IN;
}

my $dbg = 0;
foreach my $alnF (@files){ 
    my $alnIO = Bio::AlignIO->new(-file => $dir.$alnF,
				  -format=>'clustalw');
    #my $alnF =  $home.'data/IES_data/msas/alignments/cluster*_pbil-deb_pbil-deb10-local/cluster.9990.aln';
    #$home.'data/IES_data/working/genes102N.aln';
#    my $alnIO = Bio::AlignIO->new(-file => $alnF,
    #				-format=>'clustalw');
    
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
	    $id =~s/P.._scaffold.*?(P.*)T(\d+)_.*/$1G$2/;
	 #   print $id,"\t";
	    if (defined($hash{$id})){      #find if it has IES
		foreach my $ies (@{$hash{$id}}){
		    my $aaLoc = $ies->{'aaLoc'}; 
		    my $msaLoc = $alnO->column_from_residue_number($id,$aaLoc);     #find the location of IES in the alignment
		    #	  print $msaLoc,'(',$ies->{'frame'},', ',$aaLoc,', ',$ies->{'ies'},")\t";
	#	    print $msaLoc,"\t";
		    $charM{$msaLoc}{$id} = $ies->{'frame'};
		}
	    }
	 #   print "\n";
	}
	#only print if at least 1 IES is present
	next unless ((scalar(keys %charM)) > 0);
	#print "IES presence character matrix \n";
	#second pass to print
	my $charMatrixF = $alnF;
	$charMatrixF =~ s/\.aln/\.dat/;
	open OUT, '>'.$dir.$charMatrixF or die $!;
	print OUT "geneName\t";
	foreach my $state (sort {$a<=>$b} keys %charM){
	    print OUT $state,"\t";
	}
	print OUT "\n";
	foreach my $seq ($alnO->each_seq()){
	    print OUT $seq->id(),"\t";
	    my $id = $seq->id();
	    foreach my $state (sort {$a<=>$b} keys %charM){
#	    print $charM{$state},"\n";
		if (defined($charM{$state}{$id})){
		    print OUT $charM{$state}{$id},"\t";
		}else{
		    print OUT "0\t";
		}
	    }
	    print OUT "\n";
	}
	close OUT;
    }
}
