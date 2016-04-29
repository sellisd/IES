#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
# parse Gblocks output and translate results to nucleotide coordinates
my $path = $ARGV[0];

opendir(DH, $path) or die $!;

open OUT, '>', catfile($path,'gblocks.dat') or die $!;

my @files = grep { /-gb.htm/ } readdir(DH);
foreach my $file (@files){
    print $file,"\n";
    open IN, $path.$file or die $!;
    $file =~ /^cluster\.(\d+)\.aln\.fa-gb\.htm$/;
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
    my $gbf = $file;
    $gbf =~ s/-gb.htm/-gb/;
    unlink $path.$gbf;
    unlink $path.$file;
}
close OUT;
