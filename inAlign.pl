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
my $iesigP; # = '/home/dsellis/data/IES/analysis/tables/pte.iesInGenes';
my $alnPath;#  = '/home/dsellis/data/IES/analysis/msas/filtered/';
my $out;
my $usage = <<HERE;

Parses the file containing IES overlapping genes and corresponding alignment files, translates overlaps to multiple sequence alignment cooridnates

usage inAlign.pl [OPTIONS]
where OPTIONS can be:
    -iesig   :  path with .iesInGenes file
    -alnPath :  path with alignment files
    -out     :  path and file name for output
    -help|?  :  this help screen

HERE

die $usage unless (GetOptions('help|?'    => \$help,
			      'iesig=s'   => \$iesigP,
			      'alnPath=s' => \$alnPath,
			      'out=s'     => \$out
		   ));

die $usage if $help;

my $prefixesR = initF();

# read ies in gene
my %iesH;
my %geneH;
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
	my $protName = gene2prot($gene, $prefixesR);
	if(defined($iesH{$ies})){
	    push @{$iesH{$ies}{'begin'}}, $begin;
	    push @{$iesH{$ies}{'stop'}}, $stop;
	}else{
	    $iesH{$ies} = {'begin' => [$begin],
			   'stop'   => [$stop]};
	}
	if(defined($geneH{$protName})){
	    push @{$geneH{$protName}}, $ies;
	}else{
	    $geneH{$protName} = [$ies];
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

open OUT, '>', $out or die $!;
print OUT "geneFamily\tgene\tbegin\tend\ties\n";
foreach my $file (@files){
    $file =~/^cluster\.(\d+)\.nucl\.fa$/;
    my $geneFamily = $1;
    my $alnF = catfile($alnPath, $file);    
    my $alnIO = Bio::AlignIO->new(-file => $alnF,
				  -format => 'fasta');
    while(my $alnO = $alnIO->next_aln){ # open alignment and find position of ies in alignment
	my $seqNo = $alnO->num_sequences;
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $id = $seq->id();
	    # find if they have IES
	    if(defined($geneH{$id})){
		# in which location and translate location to alignment coordinates
		foreach my $ies (@{$geneH{$id}}){
		    for(my $i = 0; $i <= $#{$iesH{$ies}{'begin'}}; $i++){
			my $inGbegin = $iesH{$ies}{'begin'}[$i] + 1;
			my $inGstop  = $iesH{$ies}{'stop'}[$i];
			my $msaLocBegin = $alnO->column_from_residue_number($id, $inGbegin);
			my $msaLocStop  = $alnO->column_from_residue_number($id, $inGstop);
			print OUT $geneFamily,"\t", $id, "\t", $msaLocBegin, "\t", $msaLocStop, "\t", $ies,"\n";
		    }
		}
	    }
	}
    }
}
print "done\n";
close OUT;
