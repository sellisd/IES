#!/usr/bin/perl
use warnings;
use strict;
#read genbank file and join CDS
my $file = 'output.gnbk';
open IN, $file or die $!;
#read by feature
my $inFeature = 0;
my %geneH; #hash of gene IDs {start => A, end => B, or =>C}
#other useful information codon start, translation (this is for the whole sequence)
while(my $line = <IN>){
  chomp $line;
  my @ar = split " ", $line;
  if ($ar[0] eq 'FEATURES'){
    #entered features
    $inFeature = 1;
  }elsif($ar[0] eq 'ORIGIN'){
    #exited features
    $inFeature = 0;
  }
  if($inFeature == 1){  #inside FEATURES
    if ($ar[0] eq 'CDS'){
      my $parent = $ar[1];
      my $location = $ar[2];       #just entering a feature of interest
      my $orientation;
      my $start;
      my $end;
      if ($location =~ /complement(\d+\.\.\d+)/){
	$orientation = 'compl';
	$start = $1;
	$end = $2;
      }elsif($location =~ /\d+\.\.\d+/){
	$orientation  = '.';
	$start = $1;
	$end = $2;
      }
      if(defined($geneH{$parent})){
	if($geneH{$parent}{'orientation'} != $orientation){
	  die "ti egine re paidia!";
	}
	push @{$geneH{$parent}{'coordinates'}},$start.'..'.$end;
      }else{
	$geneH{$parent} = { 'coordinates' => \[$start.'..'.$end],
			    'orientation' => $orientation };   
      }
      #save information in hash
  }else{
    print $line,"/n";
  }
}

    #read until the parent gene is found
    #add it to the list of locations for the current parent
    #keep reading unti

close IN;

