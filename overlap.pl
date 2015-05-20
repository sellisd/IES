#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Bio::Tools::GFF;
use Getopt::Long;
my $home = '/home/dsellis/';
my $floating = 1; #by default filter floating IES
my $usage = <<HERE;

Reads
  - the output of bedtools intersect with the -wo option (which genes contain IE)
  - the .bed file with IES coordinates
  - master genbank file and produces an output file with the following columns
  - geneId
  - IESid
  - IES coordinates (begin)
  - IES coordinates (end)
  - IES length

Usage:
 overlap.pl [OPTIONS]

where options can be:
 - help|?   this help screen
 - species  species abreviation: Ppr, Pbi, Pte, Pen, Pse
 - floating if true (default) filter floating
HERE
my $help;
my $speciesAbr;

die $usage unless (GetOptions('help|?' => \$help,
			      'species=s' => \$speciesAbr,
		              'floating=i' => \$floating));
die $usage if $help;
my $dataPath = $home.'data/IES_data/';
my $iesLengthF;
my $subdir;
if($speciesAbr eq 'Ppr'){
    $subdir = 'pprimaurelia/';
}elsif($speciesAbr eq 'Pbi'){
    $subdir = 'pbiaurelia/';
    if($floating){
	$iesLengthF = 'internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.fl.gff3';
    }else{
	$iesLengthF = 'internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3';
    }
}elsif($speciesAbr eq 'Pte'){
    $subdir = 'ptetraurelia/';
    if($floating){
	$iesLengthF = 'internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3';
    }else{
	$iesLengthF = 'internal_eliminated_sequence_PGM_IES51.pt_51.gff3';
    }
}elsif($speciesAbr eq 'Pen'){
    $subdir = 'ppentaurelia/';
}elsif($speciesAbr eq 'Pse'){
    $subdir = 'psexaurelia/';
    if($floating){
	$iesLengthF = 'internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.fl.gff3';
    }else{
	$iesLengthF = 'internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3';
    }
}else{
    print 'Not known species abreviation: $speciesAbr',"\n";
    die $usage;
}

#./makeBed.pl creates bed files for IES and CDS from the IES.gnbk file
#using bedtools calculate overlap, which genes have an IES in them
# bedtools intersect -a PTET.CDS.bed -b PTET.IES.bed -wo > Pte.overlap

# parse overlap and load CDSs overlapping an IES
#$hash{$scaffold}{$CDS.id => [IES.id,...]}
my %overlapH;
$dataPath = $dataPath.$subdir;
my $overlapF = $speciesAbr.'.overlap';
my $output = $dataPath.$speciesAbr.'.IESinCDS';
open OR, $dataPath.$overlapF or die $!;
while(my $line = <OR>){
    chomp $line;
    my @ar = split " ", $line;
    my $scaffold = $ar[0];
    my $CDSiD = $ar[3];
    my $IESiD = $ar[7];
    if(defined($overlapH{$scaffold}{$CDSiD})){
	push @{$overlapH{$scaffold}{$CDSiD}},$IESiD;
    }else{
	$overlapH{$scaffold}{$CDSiD} = [$IESiD];
    }
}
close OR;

#read IES files and build
#$iesH{$iesid} = [start,end]
my $iesF = $speciesAbr.'.IES.bed';
my %iesH;

open IES, $dataPath.$iesF or die $!;
while(my $line = <IES>){
    chomp $line;
    (my $scaffold, my $start, my $end, my $iesid) = split " ", $line;
    $iesH{$scaffold}{$iesid} = [$start,$end];
}
close IES;

#Make hash with IES lengths
my %iesLengthsH;
open IESL, $dataPath.$iesLengthF or die $!;
while (my $line = <IESL>){
    chomp $line;
    my @ar = split " ", $line;
    my $scaffold = $ar[0];
    my $string = $ar[8];
    my @details = split ';', $string;
    my $id;
    my $name;
    my $length;
    foreach my $pairs (@details){
	if($pairs =~ /^ID=(.*)$/){
	    $id = $1;
	}elsif($pairs =~ /^Name=(.*)$/){
	    $name = $1;
	}elsif($pairs =~ /^sequence=(.*)$/){
	    $length = length($1);
	}
    }
    $iesLengthsH{$scaffold.$id} = $length;
}
close IESL;

#read genbank file
#and for each CDS that has an IES add up length up until the appropriate position
my $genbankF = $speciesAbr.'.IES.gnbk';
my $gnbkIn = Bio::SeqIO->new('-file' => $dataPath.$genbankF,
			     '-format' => 'genbank');
open OUT, '>'.$output or die $!;
while(my $seqO = $gnbkIn->next_seq()){
    my $scaffold = $seqO->accession_number();
    foreach my $featureO ($seqO->get_SeqFeatures()){
	if($featureO->primary_tag() eq 'CDS'){
	    my @geneName = $featureO->get_tag_values('gene');
	    if(defined($overlapH{$scaffold}{$geneName[0]})){
		my $iesInCDS = $overlapH{$scaffold}{$geneName[0]};
		#this gene has at least one IES
#usually there are few IES in each gene so checking all should not be too slow
		my $locationO = $featureO->location();
		my $complement =  $locationO->strand();
		my $index = 0;
		my @iesIndexProt; #location of IES in protein coordinates
		foreach my $locations ($locationO->each_Location()){
		    my $start = $locations->start;
		    my $end = $locations->end;
		    foreach my $iesInCDSName (@{$iesInCDS}){
			my $IESstart = @{$iesH{$scaffold}{$iesInCDSName}}[0];
			my $IESend = @{$iesH{$scaffold}{$iesInCDSName}}[1];
			if($IESstart >= $start and $IESstart <= $end){
			    #this IES is (at least partially) in this CDS
			    push @iesIndexProt, [$index+($IESstart-$start+1),$index+($IESend-$start + 1)];
			}
		    }
		    #if there is an IES between start end
		    #print gene IES name index+ (IES location - CDS start+1)
		    $index += ($end-$start+1) #length of each CDS
		}
		my $counter = 0;
		foreach my $iesIP (@iesIndexProt){
		    my $iesID = ${$iesInCDS}[$counter];
		    print OUT $geneName[0], ' ', $iesID, ' ';
		    my $gCoo; # gene coordinates
#		    my $pCoo; #protein coordinates
		    if($complement == -1){
			# if gene in - strand location is counding from the end
			$gCoo = [$index - ${$iesIP}[1] + 1, $index - ${$iesIP}[0] + 1];
#			$pCoo = $index - $iesIP; # IES has length 2
		    }elsif($complement == 1){
			$gCoo = $iesIP;
#			$pCoo = $iesIP;
		    }else{
			die $complement;
		    }
		    my $length = $iesLengthsH{$scaffold.$iesID};
#		    my $aaCoo = int(($pCoo - 1)/3) + 1; #position in amino-acid coordinates
#		    my $frame = ($pCoo - 1) % 3 + 1;
#		    print OUT $aaCoo, ' ',$frame, ' ', $length,' ',$pCoo,"\n";
		    print OUT ${$gCoo}[0], ' ',${$gCoo}[1], ' ', $length,"\n";
# nt 1 2 3 4 5 6 7
#make zero-based and add back 1 to the result
# aa 0 0 0 1 1 1 2 (a-1) / 3 + 1
# p  1 2 0 1 2 0 1 (a-1) % 3 +1
		    $counter++;
		}
	    }
	}
    }
}
close OUT;
