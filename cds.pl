#!/usr/bin/perl
use warnings;
use strict;
use Bio::Tools::GFF;

## Extract CDS locations in gene coordinates

my $gffio = Bio::Tools::GFF->new(-file => '/home/dsellis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gff3',
				 -gff_version => 3);

# my $gffio = Bio::Tools::GFF->new(-file => '/home/dsellis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gff3',
# 				 -gff_version => 3);
my %genes;
my %mRNAs;
while(my $feature = $gffio->next_feature()){
    my $type = $feature->primary_tag;
    $feature->start < $feature->end or die; # make sure start/end does not reflect strand
    if($type eq 'gene'){
	my @geneNames = $feature->get_tag_values('ID');
	die "@geneNames" if $#geneNames > 0;
	if(defined($genes{$geneNames[0]})){
	    die "@geneNames";
	}else{
	    $genes{$geneNames[0]} = [$feature->start, $feature->end, $feature->strand];
	}
    }elsif($type eq 'mRNA'){
	my @rnaNames = $feature->get_tag_values('ID');
	die "@rnaNames" if $#rnaNames > 0;
	if(defined($mRNAs{$rnaNames[0]})){
	    die "@rnaNames";
	}else{
	    my @parents = $feature->get_tag_values('Parent');
	    die @parents if $#parents > 0;
	    $mRNAs{$rnaNames[0]} = $parents[0];
	}
    }elsif($type eq 'CDS'){
	my @ids = $feature->get_tag_values('ID');
	die "@ids" if $#ids > 0;
	my @parents = $feature->get_tag_values('Parent');
	die @parents if $#parents > 0;
	my $CDSGenomicStart = $feature->start;
	my $CDSGenomicEnd = $feature->end;
#	print $parents[0],' ', , "\n";
	if(defined($genes{$mRNAs{$parents[0]}})){
	    (my $geneGenomicStart, my $geneGenomicEnd, my $strand) = @{$genes{$mRNAs{$parents[0]}}};
	    # begin < end
	    my $CDSGeneStart;
	    my $CDSGeneEnd;

	    if($strand eq '1'){
		$CDSGeneStart = $CDSGenomicStart - $geneGenomicStart + 1;
		$CDSGeneEnd   = $CDSGenomicEnd   - $geneGenomicStart + 1;
	    }elsif($strand eq '-1'){
		# maintain end > start
		$CDSGeneStart = $geneGenomicEnd - $CDSGenomicEnd   + 1;
		$CDSGeneEnd   = $geneGenomicEnd - $CDSGenomicStart + 1;
	    }else{
	     	die "Unknown strand: $strand";
	    }
	print "$ids[0] $mRNAs{$parents[0]} $CDSGenomicStart $CDSGenomicEnd $CDSGeneStart $CDSGeneEnd\n";
	}else{
	    die $parents[0];
	}
    } 
}
$gffio->close;
