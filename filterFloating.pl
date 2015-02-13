#!/usr/bin/perl
use warnings;
use strict;
#or perhaps modify iesMACseq to print for all species and then do it in R

# read F character matrices with nucleotide coordinates
# split header start.end
# calculate distances and see if two are very close (or overlapping)
# for the columns that are too close (distance 1)
# find if MAC sequence is ...TATA...
# for the IESes find if sequence is TA..TATA and TATA..TA
# if so they are floating, then merge columns and print new corrected character matrix (do not overwright though.

# read alignments/genbank file and test whether they are separated by TAs

my $minDist = 1; # if IES are less than 1pb apart it is suspicious
my $path = '/home/dsellis/data/IES_data/';
my @iesFiles = qw(
pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3
ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3
psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3
);
my $PtetMacF = $path.'ptetraurelia/Pte.ies.mac_seq';

my %IESH;
foreach my $file (@iesFiles){
    my $species = substr($file,0,3);
# parse gff3 with annotations
    open IN, $path.$file or die $!;;
    while (my $line = <IN>){
	chomp $line;
	my @ar = split "\t", $line;
	my $scaffold = $ar[0];
	my $start = $ar[3];
	my $end = $ar[4];
	my $score = $ar[5];
	my $annotation = $ar[8];
	my $id;
	my $sequence;
	my $alt_seq;
	my $mac_seq;
	my @annotations = split ';', $annotation;
	my $counter = 0;
	foreach my $annot (@annotations){
	if ($annot =~ /ID=(.*)/){
	    $id = $1;
	}elsif($annot =~ /sequence=(.*)/){
	    $sequence = $1;
	}elsif($annot =~ /alternative_IES_seq=(.*)/){
	    $alt_seq = $1;
	}elsif($annot =~ /mac_seq=(.*)/){
	    $mac_seq = $1;	    
	}
	$counter++;
	}
	unless(defined($species) and defined($scaffold) and defined($start) and defined($end) and defined($score) and defined($sequence)){
	    die $id;
	}
	my $entry = {
	    'species' => $species,
	    'scaffold' => $scaffold,
	    'start' => $start,
	    'end'   => $end,
	    'score' => $score,
	    'mac_seq' => $mac_seq,
	    'sequence' => $sequence,
	    'alt_sequence' => $alt_seq};
	if(defined($IESH{$id})){
	    die "double name! at $id\n";
	}else{
	    $IESH{$id} = $entry;
	}
    }
    close IN;
}

#load mac_seq for P. tetraurelia
open IN, $PtetMacF or die $!;
while (my $line =<IN>){
	chomp $line;
	(my $iesID, my $mac_seq) = split "\t", $line;
	$IESH{$iesID}{'mac_seq'} = $mac_seq;
}
close IN;


# foreach my $id (sort keys %IESH){
#     unless(defined($IESH{$id}{'mac_seq'})){
# 	die;
#     }
#     print $id,"\t";
#     print $IESH{$id}{'species'},"\t";
#     print $IESH{$id}{'mac_seq'},"\t";
#     print $IESH{$id}{'sequence'},"\t";
#     if (defined($IESH{$id}{'alt_sequence'})){
# 	print $IESH{$id}{'alt_sequence'},"\t";
#     }else{
# 	print "NA\t";
#     }
#     print "\n";
# }

# read character matrices
my $charM = '/home/dsellis/data/IES_data/msas/alignments/aln/cluster.2719.F.dat';
my $gBlocks = '/home/dsellis/data/IES_data/msas/alignments/aln/cluster.2719.nucl.fa.gblocks';

open GB, $gBlocks or die $!;
my @gbStart;
my @gbEnd;
while (my $line = <GB>){
    chomp $line;
    (my $start, my $end) = split " ", $line;
    push @gbStart, $start;
    push @gbEnd, $end;
}
close GB;
open CM, $charM or die $!;
my @start;
my @end;
my $line = readline(CM);
chomp $line;
my @ar = split " ", $line;
shift @ar; #remove geneName
foreach my $se (@ar){
    $se =~ /(\d+)\.(\d+)/ or die $se;
    my $start = $1;
    my $end = $2;
    push @start, $start;
    push @end, $end;
}
# do not close until we are sure we don't need to read the rest

my @within; # index of IES start whithin gblocks
my @wStart; # copies with only IES of interest
my @wEnd;
# include only IES that are not within GBlocks
for (my $i = 0; $i <= $#gbStart; $i++){
    for (my $j = 0; $j <= $#start; $j++){
	if($start[$j] >= $gbStart[$i] and $end[$j] <= $gbEnd[$i]){
	    # within
	    push @wStart, $start[$j];
	    push @wEnd, $end[$j];
	    push @within, $j;
	}
    } 
}
print "@wStart\n";
print "@wEnd\n";

# calculate distances
my %closeDist; # suspiciously close distance
for (my $i = 1; $i <= $#within; $i++){
    my $distances = $start[$within[$i]] - $end[$within[$i - 1]];
    if ($distances <= $minDist){
	$closeDist{$i} = $distances; #index of distance = value of distance
    }
}
my $distancesNo = keys %closeDist;

if ($distancesNo >= 0){
# if suspicious distances read the rest of the char matrix
# find which IESs are suspiciously close to each other
# check their MAC_sequences and IES ends
    my @upstream;
    my @downstream;
    while (my $line = <CM>){
	chomp $line;
	my @iesIDs = split " ", $line;
	my $geneName = shift @iesIDs;
	for(my $i = 0; $i <= $#iesIDs; $i++){
	    my $state = $iesIDs[$i];
	    if(defined($closeDist{$i})){ # downstream
		if($state ne '0'){  # unless absent
		    push @downstream, $state;
		}
		my $upstreamState = $iesIDs[$i-1]; # upstream
		if($upstreamState ne '0'){
		    push @upstream, $upstreamState;
		}
	    }
	}
    }
#is floating?
#if floating adjust location
    print "upstream\n";
    foreach my $upIES (@upstream){
	my $mac_seq = $IESH{$upIES}->{'mac_seq'};
	$mac_seq=~/([actg]+)([ACTG]+)[actg]+/ or die;
	my $window = length($1);
	my $ies_junction = length($2);
	my $macBorders = substr($mac_seq,$window-2, $ies_junction + 4);
	my $seq = $IESH{$upIES}->{'sequence'};
	print substr($seq,0,4), '...',substr($seq,-4,4),"\t". $macBorders,"\n";
    }
    print "downstream\n";
    foreach my $doIES (@downstream){
	my $mac_seq = $IESH{$doIES}->{'mac_seq'};
	$mac_seq=~/([actg]+)([ACTG]+)[actg]+/ or die;
	my $window = length($1);
	my $ies_junction = length($2);
	my $macBorders = substr($mac_seq,$window-2, $ies_junction + 4);
	my $seq = $IESH{$doIES}->{'sequence'};
	print substr($seq,0,4), '...',substr($seq,-4,4),"\t". $macBorders,"\n";
    }
}
close CM;

