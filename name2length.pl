#!/usr/bin/perl
use warnings;
use strict;

#read character matrix with IES names and printout character matrix with IESlength (or other properties?) as states

my $dataPath = '/home/dsellis/data/IES_data/';
my $charMatPath = $dataPath.'msas/alignments/charMat/';
# files with IES information
#read separately for ptet, pbi, pse
my @iesFiles = qw(
ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3
pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.fl.gff3
psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.fl.gff3
pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.fl.gff3
);

my %iesH;
print "reading IES information\n";
foreach my $file (@iesFiles){
    my $addTA;
    if(substr($file,0,4) eq 'pcau'){
	$addTA = 1;
    }else{
	$addTA = 0;
    }
    open IN, $dataPath.$file or die $!;
    while (my $line = <IN>){
	chomp $line;
	my @ar = split "\t", $line;
	my $scaffold = $ar[0];
	my $annotation = $ar[8];
	my $id;
	my $start = $ar[3];
	my $end = $ar[4];
	my $length;
	my $seq;
	my @annotations = split ';', $annotation;
	foreach my $annot (@annotations){
	    if ($annot =~ /ID=(.*)/){
		$id = $1;
	    }elsif($annot =~ /sequence=(.*)/){
		$seq = $1;
	    }
	}
	$length = length($seq);
	if ($addTA == 1){
	    $length+=2;
	}
	$iesH{$id} = {
	    'start'    => $start,
	    'end'      => $end,
	    'length'   => $length
		
	}
    }
    close IN;
}

# character matrices
opendir(DH, $charMatPath) or die $!;
my @charMatF = grep{/cluster\.\d+\.dat/} readdir(DH);


#read character matrix and replace
foreach my $charMF (@charMatF){
    print $charMatPath.$charMF,"\n";
    my $outputFile = $charMF;
    $outputFile =~ s/\.dat$/.L.dat/ or die;
    open IN, $charMatPath.$charMF or die $!;
    open OUT, '>'.$charMatPath.$outputFile or die $!;
    my $lineCounter = 0;
    while(my $line=<IN>){
	if($lineCounter == 0){
	    print OUT $line;
	}else{
	    chomp $line;
	    my @ar = split "\t", $line;
	    print OUT shift @ar,"\t";
	    foreach my $entry (@ar){
		my @iesEntry = split ",",$entry;
		foreach my $iesId (@iesEntry){
		    if(defined($iesH{$iesId})){
			print OUT $iesH{$iesId}{'length'},"\t";
		    }elsif($iesId eq '0'){
			print OUT '0',"\t";
		    }elsif($iesId eq 'NA'){
			print OUT 'NA',"\t";
		    }else{
			die  $iesId,"\n";
		    }
		}
	    }
	    print OUT "\n";
	}
	$lineCounter++;
    }
    close OUT;
    close IN;
}
