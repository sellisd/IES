#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use Parallel::ForkManager;

my $pm = Parallel::ForkManager->new(7);
my %reruns;

# chunks that should rerun
my $rerunF = '/home/dsellis/data/IES/analysis/allvsall/runs/2rerun.dat';

# results from cluster
my $blastoutD = '/home/dsellis/data/IES/analysis/allvsall/blastout/';

#read files that need to be rerun

open IN, $rerunF or die $!;
while(my $line = <IN>){
    $line =~ /^(\d+),.*/;
    $reruns{$1} = 1;
}
close IN;

opendir(DH, $blastoutD) or die $!;
my @nodes = readdir(DH);
close DH;
my %clusterF;
my @rerunLocally;
foreach my $node (@nodes){
    next if $node eq '.' or $node eq '..';
    my $nodePath = catfile($blastoutD, $node);
    if (-d $nodePath){
	print $nodePath,"\n";
	opendir(ND, $nodePath);
	my @files = grep {/\.dat$/} readdir(ND);
	foreach my $file (@files){
	    $file =~ /blastout\.chunk\.(\d+)\.fa\.dat/;
	    my $chunk = $1;
	    if(defined($clusterF{$file})){
		die $file," duplicated \n";
	    }else{
		if(defined($reruns{$chunk})){
		    push @rerunLocally, $chunk;
	#	    print $file," should rerun\n";
		}else{
		    $clusterF{$file} = $node
		}
	    }
	}
	close ND;
    }
}


foreach my $chunk (@rerunLocally){
    my $pid = $pm->start and next;
    my $cmdl = 'blastp -query /home/dsellis/data/IES/tempdat/fastachunks/chunk.'.$chunk.'.fa -db /home/dsellis/data/IES/analysis/protdb/allprot -outfmt 6 -out /home/dsellis/data/IES/analysis/allvsall/blastout/local/blastout.chunk.'.$chunk.'.fa.dat -seg yes';
    print $cmdl,"\n";
    system $cmdl;
    $pm->finish;
}
