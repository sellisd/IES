#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;

my %strangeExons;
my @species = qw/Pbi Pte Pse/;
foreach my $species (@species){
    my $gff3;
    if($species eq 'Pbi'){
	$gff3 = '/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_annotation_v1.gff3';
    }elsif($species eq 'Pte'){
	$gff3 = '/Users/diamantis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.4.gff3';
    }elsif($species eq 'Pse'){
	$gff3 = '/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3';
    }else{
	die;
    }
    my $gff3In = Bio::Tools::GFF->new('-file' => $gff3,
				      '-gff_version' => 3);
    
    while(my $feature = $gff3In->next_feature()){ # one line at a time
	if ($feature->primary_tag() eq 'exon'){
	    my $start = $feature->start();
	    my $end = $feature->end();
	    if ($end - $start + 1 == 102){
		my @parent = $feature->get_tag_values('Parent');
		my $parent;
		$parent = $parent[0];
		$parent =~ s/PBIGNT(\d+)/PBIGNP$1/ if $species eq 'Pbi';
		$parent =~ s/PTET.51.T(\d+)/PTET.51.P$1/ if $species eq 'Pte';
		$parent =~ s/PSEXGNT(\d+)/PSEXPNG$1/ if $species eq 'Pse';
		if(defined($strangeExons{$parent})){
		    $strangeExons{$parent}++;
		}else{
		    $strangeExons{$parent} = 1;
		}
	    }
	}
    }
}

open SLX, '/Users/diamantis/data/IES_data/working/silix.output' or die $!;
while(my $line = <SLX>){
    chomp $line;
    (my $group, my $gene) = split " ", $line;
    if(defined($strangeExons{$gene})){
	print $group,"\n";
    }
}
close SLX;

