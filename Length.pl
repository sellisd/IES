#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
my %scaffoldH;
my $scaffoldsF = $ARGV[0]; #scaffold file
my $gff3 = $ARGV[1]; #gff3 file
my $scaffoldIn = Bio::SeqIO->new('-file' => $scaffoldsF,
				 '-format' => 'fasta');
while(my $scaffoldSeq = $scaffoldIn->next_seq){
#scaffold length, #102, #101+103, #total
    $scaffoldH{$scaffoldSeq->display_id} = [$scaffoldSeq->length(),0,0,0];
}
my $gff3In = Bio::Tools::GFF->new('-file' => $gff3,
				  '-gff_version' => 3);

while(my $feature = $gff3In->next_feature()){ # one line at a time
    my $scaffold = $feature->seq_id();
    if ($feature->primary_tag() eq 'CDS'){
	if(defined($scaffoldH{$scaffold})){
	    my $length = $feature->length();
	    if($length == 102){
		${$scaffoldH{$scaffold}}[1]++;
	    }elsif($length == 101 or $length == 103){
		${$scaffoldH{$scaffold}}[2]++;
	    }
	    ${$scaffoldH{$scaffold}}[3]++;
	}
    }
}

foreach my $sc(keys %scaffoldH){
    print $sc," @{$scaffoldH{$sc}}\n";
}

