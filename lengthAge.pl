#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use File::Spec::Functions qw(catfile);
use File::HomeDir;
my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};

my $asrRun = $ARGV[0];
#my $homeD = File::HomeDir->my_home;
my $notationF =  catfile($basePath, 'analysis', 'notation.csv');

# prepare tidy table with age and length information
my $homF = catfile($basePath, 'analysis', 'iesdb', 'homIESdb'.$asrRun.'.tab');

# read homologous IES information
my %h;
open IN, $homF or die $!;
while(my $line = <IN>){
  chomp $line;
  (my $id, my $geneFamily, my $ies, my $age) = (split " ", $line)[0, 1, 9, 10];
  next if $ies eq 'NA';
  $h{$ies} = [$id, $geneFamily, $age];
}
close IN;

# read iesdb files
printab('abr', 'iesId', 'length', 'homIESId', 'geneFamily', 'age');
my $nr = getNotation($notationF);
foreach my $sp (sort keys %$nr){
  my %pab = %{$nr->{$sp}};
  #    print $pab{'abr'},"\n";
  my $f = catfile($basePath, 'analysis', 'iesdb', $pab{'abr'}.'.iesdb');
  open IN, $f or die $f;
  readline(IN);
  while(my $line = <IN>){
    chomp $line;
    (my $id, my $length) = (split " ", $line)[0, 7];
    my $iesId = $pab{'abr'}.'.'.$id;
    if(defined($h{$iesId})){
      printab($pab{'abr'}, $id, $length, @{$h{$iesId}});
    }else{
      #	    print "$iesId";
    }
  }
  close IN;
}
#iesId speciesabr length age homeiesid
