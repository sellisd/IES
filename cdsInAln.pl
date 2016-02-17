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
	foreach my $seq ($alnO->each_seq()) {     #find which genes
	    my $protId = $seq->id;
	    printab $protId;
	    my $geneId = prot2gene($protId);
	    if($geneId){
		if(defined($cdsLoc{$geneId})){
		    my @be = @{$cdsLoc{$geneId}};
		    foreach my $be (@be){
			my $l = $alnO->length;
#			use Data::Dumper;
#			my $slice = $alnO->slice(1145,1172);
# my $debug = Bio::AlignIO->new(-file => '>temp',
# 			      -format => 'fasta');
# # 			$debug->write_aln($slice)
# 			for(my $i = 1; $i <= 1140; $i++){
# 			    print $protId, ' ', $be, ' ', $alnO->column_from_residue_number($protId, $i),"\n";
# 			}
# 			die;
			
			print $alnO->column_from_residue_number($protId, $be),"\n";
			# die;
			# print $alnO->column_from_residue_number($protId, $be);
			# print ' ';
		    }
		    print "\n";
		    die;
		}
	    }
#	    my $msaLocStart = $alnO->column_from_residue_number($id,$start);     #find the 
	}
    }
}
close DH;


#     my $seqNo = $alnO->num_sequences;
# }
# }
