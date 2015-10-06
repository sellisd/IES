#!/usr/bin/perl
use warnings;
use strict;
use Math::Random::Secure qw(irand);

my $file = $ARGV[0];

my $lengthCutoff = 20;
my $identityCutoff = 95;
my $evalueCutoff = 0.001;
# blast output
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# for each IES keep the lowest e-value and if multiple the lowest bit-score
my %tophits;

open IN, $file or die $!;
while(my $line = <IN>){
    my @ar = split " ", $line;
    next if $ar[0] eq $ar[1]; # skip self hits
    #next if $ar[2] <= $identityCutoff;
    #next if $ar[3] <= $lengthCutoff;
    next if $ar[10] >= $evalueCutoff;
    my $currentMatch = {'sseqid' => $ar[1],
			'pident' => $ar[2],
			'length' => $ar[3],
			'mismatch' => $ar[4],
			'gapopen' => $ar[5],
			'qstart' => $ar[6],
			'qend' => $ar[7],
			'sstart' => $ar[8],
			'send' => $ar[9],
			'evalue' => $ar[10],
			'bitscore' => $ar[11],
			'out of' => 0};
    if(defined($tophits{$ar[0]})){ # is this IES previously seen?
	if($tophits{$ar[0]}->{'evalue'} < $ar[10]){ # if evalue is less
	    $tophits{$ar[0]} = $currentMatch;       #     replace
	}elsif($tophits{$ar[0]}->{'evalue'} == $ar[10]){  # if equal
	    if($tophits{$ar[0]}->{'bitscore'} > $ar[11]){ # compare bitscore
		$tophits{$ar[0]} = $currentMatch;
	    }elsif($tophits{$ar[0]}->{'bitscore'} == $ar[11]){
#if same e-value and bitscore
#pick one randomly
#		print $ar[0], $ar[1];
#		die;
		my $int = irand(2);
		if($int == 0){
		    $currentMatch->{'out of'} = $tophits{$ar[0]}->{'out of'};
		    $tophits{$ar[0]} = $currentMatch;
		}
		$tophits{$ar[0]}->{'out of'}++;
	    }
	}
    }else{
	$tophits{$ar[0]} = $currentMatch;
    }
}
close IN;

foreach my $qies (sort keys %tophits){
    print &iesName2species($qies), ' ';
    print &iesName2species($tophits{$qies}->{'sseqid'}),' ';
    print $tophits{$qies}->{'pident'},' ';
    print $tophits{$qies}->{'length'},' ';
    print $tophits{$qies}->{'mismatch'},' ';
    print $tophits{$qies}->{'evalue'},' ';
    print $tophits{$qies}->{'bitscore'},' ';
    print $tophits{$qies}->{'out of'},' ';
    print "\n";
}


sub iesName2species{
#1 P. biaurelia
#2 P. tetraurelia
#3 P. sexaurelia
#4 P. caudatum
#-1 error/unknown
    my $iesName = shift @_;
    my $sp = substr($iesName,0,11);
    if($sp eq 'IESMIC.PBIA'){
	return 1;
    }elsif($sp eq 'IESPGM.PTET'){
	return 2;
    }elsif($sp eq 'IESMIC.PSEX'){
	return 3;
    }elsif(substr($sp,0,4) eq 'MICA'){
	return 4;
    }else{
	return -1;
    }
}
