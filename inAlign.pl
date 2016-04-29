#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Bio::AlignIO;
use File::Spec::Functions qw(catfile);
use Getopt::Long;

# option variables
my $help;
my $iesigP;
my $alnPath;
my $out;
my $usage = <<HERE;

Parses the file containing IES overlapping genes and corresponding alignment files, translates overlaps to multiple sequence alignment cooridnates

usage inAlign.pl [OPTIONS]
where OPTIONS can be:
    -iesig   :  path with .iesInGenes file
    -alnPath :  path with alignment files
    -out     :  path and file name for output (.tab and .be will appended for each file)
    -help|?  :  this help screen

HERE

die $usage unless (GetOptions('help|?'    => \$help,
			      'iesig=s'   => \$iesigP,
			      'alnPath=s' => \$alnPath,
			      'out=s'     => \$out
		   ));

die $usage if $help;

my $prefixesR = initF();

my $outTab = $out.'.tab';
my $outBed = $out.'.be';

# read ies in gene
my %iesH;
my %geneH; # hash with IESs in genes
print "reading ies in gene data ... ";
opendir IDH, $iesigP or die "$! $iesigP";
my @iesigF = grep {/^.*\.iesInGenes$/} readdir(IDH);
close IDH;

foreach my $iesigF (@iesigF){
    print $iesigF,' ';
    $iesigF =~ /^(\w\w\w)\.iesInGenes$/;
    my $abr = $1;
    open IG, catfile($iesigP, $iesigF) or die $!;
    while(my $line = <IG>){
	(my $gene, my $begin, my $stop, my $ies) = split " ", $line;
	$ies = $abr.'.'.$ies;
	$begin++; # from 0 based indexing to 1-based
	my $protName = gene2prot($gene, $prefixesR);
	die "$line" if $begin == 0;
	if(defined($iesH{$ies})){
	    push @{$iesH{$ies}{'begin'}}, $begin;
	    push @{$iesH{$ies}{'stop'}}, $stop;
	}else{
	    $iesH{$ies} = {'begin' => [$begin],
			   'stop'   => [$stop]};
	}
	if(defined($geneH{$protName})){
	    $geneH{$protName}{$ies} = 1;
	}else{
	    $geneH{$protName} = {$ies =>1};
	}
    }
    close IG;
}
print "done\n";

#once all are read go trhough alignments one by one
print "reading alignments ...";
opendir DH, $alnPath or die "$! $alnPath";
my @files = grep {/\.nucl\.fa$/} readdir(DH);
close DH;

open OUT, '>', $outTab or die $!;
open BE, '>', $outBed or die $!;
print OUT "geneFamily\tgene\tbegin\tend\ties\n";
my $fileCounter = 0;
foreach my $file (@files){
    print $fileCounter, $#files, "\r";
    $file =~/^cluster\.(\d+)\.nucl\.fa$/;
    my $geneFamily = $1;
    my $alnF = catfile($alnPath, $file);    
    my $alnIO = Bio::AlignIO->new(-file => $alnF,
				  -format => 'fasta');
    while(my $alnO = $alnIO->next_aln){ # open alignment and find position of ies in alignment
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $id = $seq->id();
	    # find if they have IES
	    if(defined($geneH{$id})){
		my %iesRanges;
		# in which location and translate location to alignment coordinates
		foreach my $ies (keys %{$geneH{$id}}){
		    for(my $i = 0; $i <= $#{$iesH{$ies}{'begin'}}; $i++){
			my $inGbegin = $iesH{$ies}{'begin'}[$i];
			my $inGstop  = $iesH{$ies}{'stop'}[$i];
			if ($inGbegin == 0){
			    die;
			}
			my $msaLocBegin = $alnO->column_from_residue_number($id, $inGbegin);
			my $msaLocStop  = $alnO->column_from_residue_number($id, $inGstop);
			print OUT $geneFamily,"\t", $id, "\t", $msaLocBegin, "\t", $msaLocStop, "\t", $ies,"\n";
			if(defined($iesRanges{$ies})){
			    if($msaLocBegin < $iesRanges{$ies}[0]){
				$iesRanges{$ies}[0] = $msaLocBegin;
			    }
			    if($msaLocStop > $iesRanges{$ies}[1]){
				$iesRanges{$ies}[1] = $msaLocStop;
			    }
			}else{
			    $iesRanges{$ies} = [$msaLocBegin, $msaLocStop];
			}
		    }
		    my $beStart = $iesRanges{$ies}[0];
		    $beStart--; # from 0-indexed to 1-indexed
		    my $beEnd = $iesRanges{$ies}[1];
		    print BE $geneFamily,"\t", $beStart, "\t", $beEnd, "\t", $ies, "\t", $id, "\n";
		}
	    }
	}
    }
    $fileCounter++;
}
print "done\n";
close OUT;
close BE;
