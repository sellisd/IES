#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
my %hash;
my @files;
my $dir = $ARGV[0];
opendir(DH, $dir) or die $!;
@files = grep { /aln/ } readdir(DH);
#print @files," \n";die;
my @iesF = ('/Users/dsellis/data/IES_data/pbiaurelia/Pbi.IESinCDS',
	    '/Users/dsellis/data/IES_data/ptetraurelia/Pte.IESinCDS',
	    '/Users/dsellis/data/IES_data/psexaurelia/Pse.IESinCDS');
foreach my $file (@iesF){
  open IN, $file or die $!;
  while (my $line = <IN>){
    chomp $line;
    (my $gene, my $ies, my $aaLoc, my $frame) = split " ",$line;
    $gene =~ /(.*)G(\d+)/;
    my $speciesAbr = $1;
    my $number = $2;
    $gene = $speciesAbr.'P'.$number; #to match the proteins
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
 use Data::Dumper;
 #print Dumper %hash;
 #die;
my $dbg = 0;
foreach my $alnF (@files){ 
  my $alnIO = Bio::AlignIO->new(-file => $dir.$alnF,
				-format=>'clustalw');
  
  while(my $alnO = $alnIO->next_aln){
    my $seqNo = $alnO->num_sequences;
    print $alnF,' ',$seqNo,"\t";
    print $alnO->percentage_identity(),"\n";
    next if ($seqNo < 3);
    #for each sequence in alignment
    foreach my $seq ($alnO->each_seq()) {     #find which genes
      #do something with $seq
      my $id = $seq->id();
      print $id,"\t";
      if (defined($hash{$id})){      #find if it has IES
	foreach my $ies (@{$hash{$id}}){
	  my $aaLoc = $ies->{'aaLoc'}; 
	  my $msaLoc = $alnO->column_from_residue_number($id,$aaLoc);     #find the location of IES in the alignment
	  print $msaLoc,'(',$ies->{'frame'},', ',$aaLoc,', ',$ies->{'ies'},")\t";
	}
      }
      print "\n"
    }

    #make character matrix and print
  }
  if ($dbg>10){die;}
$dbg++;
}
  #read IESinCDS file build hash
#read alignment file(s)
#for each sequence find if there are IES and which are the locations (and frames)
#translate location in alignment location
#printout



##auto tha einai ena arxeio xehoristo
#1. overlap.pl
# read IES locations
# read CDS locations
# for each CDS find if it contains an IES and where in the protein coordinates it would be
# print protein, IES location lists as genebank? format



#kai edo ena deutero gia auti ti douleia
#read a protein MSA and find if proteins contain IES and at which locatiosn
#find correpsonding location in MSA get_seq_by(pos) Bio::SimpleAlign
#print MSA with some annotation for IES

#2. if an IES is not in a gene, what is the distance from the nearest gene, how does this compare to random points?
