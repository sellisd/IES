#!/usr/bin/perl
use warnings;
use strict;

#read genbank file and join CDS
#print lines before features
#then parse features and join CDS the rest push in an array

my $file = 'output.gnbk';
open IN, $file or die $!;
#read by feature
my $inFeature = 0;
my %geneH; #hash of gene IDs {start => A, end => B, or =>C}
#other useful information codon start, translation (this is for the whole sequence)
my $curGene;
my @coordinates;
my $curOrientation;
while(my $line = <IN>){
  chomp $line;
  if($line eq '//'){
    #reset
    $inFeature = 0;
    $curGene = undef;
    @coordinates = undef;
    $curOrientation = undef;
  }
  my @ar = split " ", $line;
  if ($ar[0] eq 'FEATURES'){
    #entered features
    $inFeature = 1;
    print $line,"\n";
    next;
  }elsif($ar[0] eq 'ORIGIN'){
    #exited features
    $inFeature = 0;
    print $line,"\n";
    next;
  }
  if($inFeature == 1){  #inside FEATURES
####
parse features









    if($ar[0] eq 'gene'){
      print $line,"\n"; # print gene start..end
      $line = <IN>;     #
      print $line; # print      /gene="XXX"
      while(1){
	#keep reading CDS
	$line = <IN>;
#	chomp $line;
	my @innerAr = split " ", $line;
	if ($innerAr[0] ne 'CDS'){
	  seek(IN,-length($line),1); #go back one line;
	  last;
	}
	my $parent = $innerAr[1];
	my $location = $innerAr[2];       #just entering a feature of interest
	my $orientation;
	my $start;
	my $end;
	if ($location =~ /complement\((\d+)\.\.(\d+)\)/){
	  $orientation = 'compl';
	  $start = $1;
	  $end = $2;
	}elsif($location =~ /(\d+)\.\.(\d+)/){
	  $orientation  = '.';
	  $start = $1;
	  $end = $2;
	}
	push @coordinates, ($start, $end);
	$curOrientation =  $orientation;	
      }
      #make a strign with output and then cut it to fit in 80char width
       print "     CDS             ";
      # print "complement(" if ($curOrientation eq 'compl');
      # print "join(" if ($#coordinates > 1);
      # for(my $i = 0; $i <= ($#coordinates/2)+1; $i+=2){
      # 	print $coordinates[$i],'..',$coordinates[$i+1];
      # 	print ',' unless ($i+1 == $#coordinates); #last element
      # }
      # print ")" if ($#coordinates > 1);
      # print ")" if ($curOrientation eq 'compl');
       print "\n" ;#joined coordinates
      
      @coordinates = ();
    }else{
      print $line,"\n";		      
    }
  }else{
    print $line,"\n";		
  }
}
#use Data::Dumper; print Dumper %geneH;
#read until the parent gene is found
#add it to the list of locations for the current parent
#keep reading unti

close IN;
  
