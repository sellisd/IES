#!/usr/bin/perl
use warnings;
use strict;

my $path = $ARGV[0];

opendir(DH, $path) or die $!;
my @files = grep { /-gb.htm/ } readdir(DH);
foreach my $file (@files){
    open IN, $path.$file or die $!;
    my $output = $file;
    $output =~ s/-gb\.htm/.gblocks/ or die $!; # make sure it changes the name
    open OUT, '>'.$path.$output or die $!;
    while(my $line = <IN>){
	if(substr($line,0,7) eq 'Flanks:'){
	    chomp $line;
	    my @blocks = split /[\[\]]+/, $line;
	    shift @blocks; #remove flanks header
	    foreach my $block (@blocks){
		next if $block =~ /^\s+$/;
		(my $start, my $end) =  split " ", $block;
		if(defined($start) and defined($end)){
		    print OUT $start,"\t",$end,"\n";
		}else{
		    die "$file: $block\n $line";
		}
	    }
	}
    }
    close OUT;
    close IN;
    unlink $path.$file.'-gb';
    unlink $path.$file.'-gb.htm';
}