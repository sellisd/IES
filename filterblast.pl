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
my %worsehits;
open IN, $file or die $!;
while(my $line = <IN>){
    my @ar = split " ", $line;
    #next if $ar[2] <= $identityCutoff;
    #next if $ar[3] <= $lengthCutoff;
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
    if($ar[0] eq $ar[1]){
#pick the worse self-hits and print separately
	if(defined($worsehits{$ar[0]})){
#	    keep the onw with the worse evalue
	    if($worsehits{$ar[0]}->{'evalue'} < $ar[10]){
		$worsehits{$ar[0]} = $currentMatch;
	    }
	}else{
	    $worsehits{$ar[0]} = $currentMatch;
	}
	next;  # skip self hits
    }
    next if $ar[10] >= $evalueCutoff; # keep best hits with high e-value only
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

# foreach my $qies (sort keys %worsehits){
#     print $qies,' ';
#     print $worsehits{$qies}->{'sseqid'},' ';
#     print $worsehits{$qies}->{'pident'},' ';
#     print $worsehits{$qies}->{'length'},' ';
#     print $worsehits{$qies}->{'mismatch'},' ';
#     print $worsehits{$qies}->{'evalue'},' ';
#     print $worsehits{$qies}->{'bitscore'},' ';
#     print $worsehits{$qies}->{'out of'},' ';
#     print "\n";
# }

foreach my $qies (sort keys %tophits){
    print &iesName2species($qies), ' ';                        # 1 query IES
    print &iesName2species($tophits{$qies}->{'sseqid'}),' ';   # 2 top hit IES
    print $tophits{$qies}->{'pident'},' ';                     # 3 pident
    print $tophits{$qies}->{'length'},' ';                     # 4 length
    print $tophits{$qies}->{'mismatch'},' ';                   # 5 mismatch
    print $tophits{$qies}->{'evalue'},' ';                     # 6 e-value
    print $tophits{$qies}->{'bitscore'},' ';                   # 7 bit score
    print $tophits{$qies}->{'out of'},' ';                     # 8 No of equally best hits
    if(defined($worsehits{$qies})){
	# The self-hit for this IES with the worse e-value
	print $worsehits{$qies}->{'evalue'},' ';               # 9 e-value of worse self-hit
    }else{
	die;
# There is no self-hit for this IES???
    }
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
