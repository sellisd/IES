#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;

# find in genes with multiple exons with length 102pb are

#   find for each gene the number of strange exons
my %strangeExons;
my @species = qw/Pbi Pte Pse/;
foreach my $species (@species){
    my $gff3;
    if($species eq 'Pbi'){
	$gff3 = '/Users/diamantis/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gff3';
    }elsif($species eq 'Pte'){
	$gff3 = '/Users/diamantis/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gff3';
    }elsif($species eq 'Pse'){
	$gff3 = '/Users/diamantis/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.gff3';
    }else{
	die;
    }
    my $gff3In = Bio::Tools::GFF->new('-file' => $gff3,
				      '-gff_version' => 3);
    
    while(my $feature = $gff3In->next_feature()){ # one line at a time
	if ($feature->primary_tag() eq 'CDS'){
	    my $start = $feature->start();
	    my $end = $feature->end();
	    if ($end - $start + 1 == 102){
		my @parent = $feature->get_tag_values('Parent');
		my $parent;
		$parent = $parent[0];
#at v2.0 these changed		# $parent =~ s/PBI.V1_4.1.G(\d+)/PBIGNP$1/ if $species eq 'Pbi';
		# $parent =~ s/PTET.51.T(\d+)/PTET.51.P$1/ if $species eq 'Pte';
		# $parent =~ s/PSEXGNT(\d+)/PSEXPNG$1/ if $species eq 'Pse';
		if(defined($strangeExons{$parent})){
		    $strangeExons{$parent}++;
		}else{
		    $strangeExons{$parent} = 1;
		}
	    }
	}
    }
}
#   for each group find how many genes with strange exons and how many total
#loop through group gene allele
#load genes in memory
#load proteins in memory
foreach my $gene (keys %strangeExons){
    print $gene,' ',$strangeExons{$gene},"\n";
    if($strangeExon{$gene}>3){
	#print N
#print P
    }
}

# my %countExons;
# my %geneGroups;
# open SLX, '/Users/diamantis/data/IES_data/working/silix.output' or die $!;
# while(my $line = <SLX>){
#     chomp $line;
#     (my $group, my $gene) = split " ", $line;
#     if(defined($geneGroups{$group})){
# 	$geneGroups{$group}++;
#     }else{
# 	$geneGroups{$group}=1;
#     }
#     if(defined($strangeExons{$gene})){
# 	if(defined($countExons{$group})){
# 	    $countExons{$group}++;
# 	}else{
# 	    $countExons{$group}=1;
# 	}
#     }
# }
# close SLX;

# foreach my $group (sort {$a<=>$b} keys %geneGroups){
#     print $group,"\t";
#     if(defined($countExons{$group})){
# 	print $countExons{$group};
#     }else{
# 	print 0;
#     }
#     print "\t";
#     print $geneGroups{$group},"\n";
# }

