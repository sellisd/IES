#!/usr/bin/perl
use warnings;
use strict;
#use Parallel::ForkManager;
use File::Spec::Functions qw(catfile);

# run Gblocks on the protein alignments and parse output

my $path = $ARGV[0];

#my $pm = Parallel::ForkManager->new(1);

open OUT, '>', catfile($path,'gblocks.dat') or die $!;

opendir(DH, $path) or die $!;

my @files = grep { /cluster\..*\.aln\.fa$/ } readdir(DH);
foreach my $file (@files){
#    my $pid = $pm->start and next;
    print $file,"\n";
    my $cmdl = 'Gblocks '.$path.$file.' -t="protein"';
    my $htmlfile = $file.'-gb.htm';
    print $cmdl,"\n";
    system $cmdl;
    my $pnf = catfile($path, $htmlfile);
    open IN, $pnf or die "$! $pnf";
    $file =~ /^cluster\.(\d+)\.aln\.fa$/;
    my $cluster = $1;
    while(my $line = <IN>){
	if(substr($line,0,7) eq 'Flanks:'){
	    chomp $line;
	    my @blocks = split /[\[\]]+/, $line;
	    shift @blocks; #remove flanks header
	    foreach my $block (@blocks){
		next if $block =~ /^\s+$/;
		(my $start, my $end) =  split " ", $block;
                #transform start and end to nucleotide coordinates
		$start = $start*3 -2;
		$end = $end*3;
		if(defined($start) and defined($end)){
		    print OUT $cluster,"\t",$start,"\t",$end,"\n";
		}else{
		    die "$file: $block\n $line";
		}
	    }
	}
    }
    close IN;
#    $pm->finish;
}
#$pm->wait_all_children;
close OUT;
