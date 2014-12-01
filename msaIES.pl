#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
my %hash;
my @files;
#my $dir = $ARGV[0];
#opendir(DH, $dir) or die $!;
#@files = grep { /aln/ } readdir(DH);

#print @files," \n";die;
my @iesF = ('/Users/diamantis/data/IES_data/pbiaurelia/Pbi.IESinCDS',
	    '/Users/diamantis/data/IES_data/ptetraurelia/Pte.IESinCDS',
	    '/Users/diamantis/data/IES_data/psexaurelia/Pse.IESinCDS');
foreach my $file (@iesF){
  open IN, $file or die $!;
  while (my $line = <IN>){
    chomp $line;
    (my $gene, my $ies, my $aaLoc, my $frame) = split " ",$line;
    $gene =~ /(.*)G(\d+)/;
    my $speciesAbr = $1;
    my $number = $2;
    $gene = $speciesAbr.'G'.$number; # replace with P for proteins 
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
#use Data::Dumper; print Dumper %hash;die;
my $dbg = 0;
#foreach my $alnF (@files){ 
#  my $alnIO = Bio::AlignIO->new(-file => $dir.$alnF,
#				-format=>'clustalw');
my $alnF =  '/Users/diamantis/data/IES_data/working/genes102N.aln';
my $alnIO = Bio::AlignIO->new(-file => $alnF,
				-format=>'clustalw');

 
  my %frameH;
  my @characterL;
  my @characterN;
  while(my $alnO = $alnIO->next_aln){
    my $seqNo = $alnO->num_sequences;
    print $alnF,' ',$seqNo,"\t";
    print $alnO->percentage_identity(),"\n";
    next if ($seqNo < 3);
    #for each sequence in alignment
    foreach my $seq ($alnO->each_seq()) {     #find which genes
      my $id = $seq->id();
      $id =~s/P.._scaffold.*?(P.*)T(\d+)_.*/$1G$2/;
      print $id,"\t";
      if (defined($hash{$id})){      #find if it has IES
	foreach my $ies (@{$hash{$id}}){
	  my $aaLoc = $ies->{'aaLoc'}; 
	  my $msaLoc = $alnO->column_from_residue_number($id,$aaLoc);     #find the location of IES in the alignment
#	  print $msaLoc,'(',$ies->{'frame'},', ',$aaLoc,', ',$ies->{'ies'},")\t";
	  print $msaLoc,"\t";
	  if(defined($frameH{$id})){
	      push @{$frameH{$id}}, $ies->{'frame'};
	  }else{
	      $frameH{$id}=[$ies->{'frame'}];
	  }
	  push @characterL,$msaLoc;
	  push @characterN,$ies->{'ies'};
	}
      }
      print "\n";
    }

    # use Data::Dumper;
    # print Dumper %frameH;
#frame
#length
    #make character matrix and print
#     print "@characterN\n";
#     foreach my $seq ($alnO->each_seq()) {     #find which genes
# 	my $id = $seq->id();
# 	print $id,' ';
# 	foreach my $charL (@characterL){
# 	    if(defined($frameH{$id})){
		
# 	    }else{

# 	    }
# 	}
#     }
#     for each sequence $frame{$sequence} = [iesloc]
# 	push @characters, $iesloc;
# foreach sequence
# foreach location
#   print if present

  }
  if ($dbg>10){die;}
$dbg++;
#}
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
