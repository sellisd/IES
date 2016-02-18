#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use functions;
use Bio::AlignIO;
use File::Spec::Functions qw(catfile);

## find location of CDS in alignment coordinates

my $alnDir = $ARGV[0];

# build a hash of cds coordinates in gene coordinates
my @cdsFiles = qw#
/home/dsellis/data/IES_data/ptetraurelia/pte.cds.prot
/home/dsellis/data/IES_data/pbiaurelia/pbi.cds.prot
#;

my %cds; # make sure cdsIds are unique
my %cdsLoc; # keep cds locations
foreach my $cdsFile (@cdsFiles){
    open IN, $cdsFile or die $!;
    my $header = readline(IN);
    while(my $line = <IN>){
	chomp $line;
	(my $cdsid, my $geneid, my $geneStart, my $geneEnd, my $protStart, my $protEnd) = split "\t", $line;
	if(defined($cds{$cdsid})){
	    die("two CDSs with the same name: $cdsid");
	}
	if(defined($cdsLoc{$geneid})){
	    push @{$cdsLoc{$geneid}}, ($protStart, $protEnd);
	}else{
	    $cdsLoc{$geneid} = [$protStart, $protEnd];
	}
    }
    close IN;
}
opendir(DH, $alnDir) or die $!;

# calculate CDS proteinStart and CDS proteinEnd from CDSgeneStart and CDSgeneEnd and gene length
my @geneFamiliesFiles = grep { /^.*\.nucl\.fa$/ } readdir(DH);

foreach my $file (@geneFamiliesFiles){
    my $alnIO = Bio::AlignIO->new(-file => catfile($alnDir, $file),
				  -format=>'fasta'); # for the nucleotide alignments
    while(my $alnO = $alnIO->next_aln){
	printab $file;
	my $outfile = $file;
	$outfile =~ s/\.nucl\.fa$/.cds.boundaries/ or die;
	$outfile = catfile($alnDir, $outfile);
	open OUT, '>'.$outfile or die $!;
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $protId = $seq->id;
	    my $geneId = prot2gene($protId);
	    if($geneId){
		if(defined($cdsLoc{$geneId})){
		    my @be = @{$cdsLoc{$geneId}};
		    print $file,"\n";
		    for(my $c = 0; $c <= $#be - 1; $c+=2){
			print OUT $geneId,"\t";
			print OUT $alnO->column_from_residue_number($protId, $be[$c]),"\t";
			print OUT $alnO->column_from_residue_number($protId, $be[$c + 1]),"\t";
			print OUT "\n";
		    }
		}
	    }
	}
	close OUT;
    }
}
close DH;
