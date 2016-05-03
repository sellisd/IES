#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use Parallel::ForkManager;

my $pm        = Parallel::ForkManager->new(7);
my $blastoutD = '/home/dsellis/data/IES/analysis/allvsall/blastout/';
my $fastaD    = '/home/dsellis/data/IES/tempdat/fastachunks/';

opendir DH, $fastaD or die $!;
my @chunkF = grep {/^chunk\.\d+\.fa$/} readdir(DH);
close DH;

foreach my $chunk (@chunkF){
    my $pid = $pm->start and next;
    my $cmdl = 'blastp -query '.catfile($fastaD, $chunk).' -db /home/dsellis/data/IES/analysis/protdb/allprot -outfmt 6 -out /home/dsellis/data/IES/analysis/allvsall/blastout/local/blastout.'.$chunk.'.fa.dat -seg yes';
    print $cmdl,"\n";
    system $cmdl;
    $pm->finish;
}
