#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

# make a hash to translate phyldog to rb node notation
open DICT,'/home/dsellis/data/IES/analysis/tables/nodeDictionary.dat' or die $!;
my %phyldog2rb;
my %geneFamilies; # gene families in dictionary

readline(DICT); #header
while(my $line = <DICT>){
    chomp $line;
    (my $cluster, my $r, my $phyldog, my $rb) = split " ", $line;
    $phyldog2rb{$cluster.'.'.$phyldog} = $rb;
    $geneFamilies{$cluster} = 1;
}
close DICT;

# make a hash to translate from rbnode id to presence probability per ies column
my %asr;
open ASR, '/home/dsellis/data/IES/analysis/tables/avNodeProb.dat' or die $!;
while(my $line = <ASR>){
    chomp $line;
    (my $cluster, my $rb, my $iesColumn, my $presence) = split " ", $line;
    if(defined($asr{$cluster.'.'.$rb})){
	$asr{$cluster.'.'.$rb}{$iesColumn} = $presence;
    }else{
	$asr{$cluster.'.'.$rb} = {$iesColumn => $presence};
    }
}

open NP, '/home/dsellis/data/IES/analysis/tables/nodePaths.dat' or die $!;

my $lineCounter = 0;
print join("\t",(qw/cluster iesColumn from to Panc gain loss/)), "\n";
while(my $line = <NP>){
    if ($lineCounter == 0){
    }else{
	chomp $line;
	(my $cluster, my $from, my $to, my $path) = split " ", $line;
	unless($geneFamilies{$cluster}){ # if we have no data for a gene family skip it
	    #print "skipping $cluster\n";
	    next;
	}
	my @pathP = split ",", $path;
	my @pathRB;
	foreach my $nodeP (@pathP){
	    push @pathRB, $phyldog2rb{$cluster.'.'.$nodeP};
	}
	my $iesColumns =  keys $asr{$cluster.'.'.$pathRB[0]};
	for(my $i = 0; $i < $iesColumns; $i++){
	    my @probs;
	    foreach my $node (@pathRB){
		push @probs, $asr{$cluster.'.'.$node}{$i};
	    }
	    (my $gain, my $loss) = &gainLoss(@probs);
	    print join("\t",($cluster, $i,$from,$to,$probs[0],$gain,$loss)),"\n";
	}
    }
    $lineCounter++;
}

close NP;
close ASR;

sub gainLoss{
    my @l = @_;
    my $gain = 0;
    my $loss = 0;
    my $current = shift @l;
    while(@l){
	my $new = shift @l;
	my $diff = $new - $current;
	if($diff > 0){
	    $gain += $diff;
	}elsif($diff < 0){
	    $loss += -$diff;
	}
	$current = $new;
    }
    return($gain, $loss);
}
