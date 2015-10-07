#!/usr/bin/perl
use warnings;
use strict;
use Math::Random::Secure qw(irand);

my $inputFile =  '/home/dsellis/data/IES_data/iesseq/iesseq.blastout';
my $outputFile = '/home/dsellis/data/IES_data/iesseq/iesseq.filt.blastout';

my $evalueCutoff = 0.001;
# blast output
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# for each IES keep the lowest e-value and if multiple the lowest bit-score
my %tophits;
my %selfhits;
open IN, $inputFile or die $!;
while(my $line = <IN>){
    my @ar = split " ", $line;
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
	 $selfhits{$ar[0]} = 1;
 	next;  # skip self hits
     }
#    next if $ar[10] >= $evalueCutoff; # keep best hits with high e-value only
    if(defined($tophits{$ar[0]})){ # is this IES previously seen?
	if($ar[10] < $tophits{$ar[0]}->{'evalue'}){ # if evalue is less
	    $tophits{$ar[0]} = $currentMatch;       #     replace
	}elsif($ar[10] == $tophits{$ar[0]}->{'evalue'}){  # if equal
	    if($ar[11] > $tophits{$ar[0]}->{'bitscore'}){ # compare bitscore
		$tophits{$ar[0]} = $currentMatch;
	    }elsif($ar[11] == $tophits{$ar[0]}->{'bitscore'}){
                #if same e-value and bitscore
                #pick one randomly
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

foreach my $qies (sort keys %selfhits){
    if(!defined($tophits{$qies})){ # is there an IES with only self-hits?
	print $qies,"\n";
    }
}

open OUT, '>'.$outputFile or die $!;
foreach my $qies (sort keys %tophits){
    print OUT &iesName2species($qies), ' ';                        # 1 query IES
    print OUT &iesName2species($tophits{$qies}->{'sseqid'}),' ';   # 2 top hit IES
    print OUT $tophits{$qies}->{'pident'},' ';                     # 3 pident
    print OUT $tophits{$qies}->{'length'},' ';                     # 4 length
    print OUT $tophits{$qies}->{'mismatch'},' ';                   # 5 mismatch
    print OUT $tophits{$qies}->{'evalue'},' ';                     # 6 e-value
    print OUT $tophits{$qies}->{'bitscore'},' ';                   # 7 bit score
    print OUT $tophits{$qies}->{'out of'},' ';                     # 8 No of equally best hits
    print OUT $qies,' ';                                           # 9 IES id
    print OUT "\n";
}
close OUT;

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
