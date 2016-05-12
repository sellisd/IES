#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use Parallel::ForkManager;

# gather all blast runs from cluster and local and combine results for SiliX

# load local copies
# load cluster nodes
# copy cluster nodes that were not run locally and check if they are unique
my $baseP = '/home/dsellis/data/IES/analysis/allvsall/blastout';
my $localP = 'local';
my @clusterP =qw/*.dat_pbil-deb_pbil-deb11-local
*.dat_pbil-deb_pbil-deb12-local
*.dat_pbil-deb_pbil-deb13-local
*.dat_pbil-deb_pbil-deb6-local
*.dat_pbil-deb_pbil-deb8-local
/;

my $out = '/home/dsellis/data/IES/analysis/allvsall/blastout/4silix.dat';
my @in;

# check local files
my %locals;
opendir(DH, catfile($baseP, $localP)) or die $!;
my @lF = grep {/^blastout\.chunk\.(\d+)\.fa\.dat$/} readdir(DH);
foreach my $f (@lF){
    if(defined($locals{$f})){
	die "double local file $f";
    }else{
	$locals{$f} = 1;
	push @in, catfile($baseP, $localP, $f);
    }
}
close DH;

# check files on cluster
my %clusterF;
foreach my $nodeP (@clusterP){
    opendir(DH, catfile($baseP, $nodeP)) or die $!;
    my @nodeF = grep {/^blastout\.chunk\.(\d+)\.fa\.dat$/} readdir(DH);
    foreach my $f (@nodeF){
	if(defined($clusterF{$f})){
	    die "double node file $f";
	}else{
	    $clusterF{$f} = 1;
	}
	if(defined($locals{$f})){
	    # if also run locally use the local version
	}else{
	    push @in, catfile($baseP, $nodeP, $f);
	}
    }
    close DH;
}
close OUT;

# merge blast output
open OUT, '>', $out or die $!;
foreach my $f (@in){
    print $f,"\n";
    open IN, $f or die $!;
    my @lines = <IN>;
    close IN;
    print OUT @lines;
}
close OUT;
