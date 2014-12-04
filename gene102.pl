#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;

#find genes in a silix cluster and extract them to a fasta file

#use Parallel::ForkManager;
#my $cores = 4;
#my $pm = Parallel::ForkManager->new($cores);

# find in genes with multiple exons with length 102pb are

#   find for each gene the number of strange exons
my $home = '/home/dsellis/';
my %strangeExons;
my @species = qw/Pbi Pte Pse/;
my $proteinF;
my $geneF;
print "# reading gff files\n";
# for(my $p = 1; $p <4; $p++){
#     my $pid = $pm->start and next;
#     $pm->finish;
# }
# $pm->wait_all_children;

my %nseq;
my %pseq;

foreach my $species (@species){
    my $gff3;
    if($species eq 'Pbi'){
	$gff3 = $home.'data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gff3';
	$proteinF = $home.'data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.protein.fa';
	$geneF = $home.'data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gene.fa';
    }elsif($species eq 'Pte'){
	$gff3 = $home.'data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gff3';
	$proteinF = $home.'data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.protein.fa';
	$geneF = $home.'data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gene.fa';
    }elsif($species eq 'Pse'){
	$gff3 = $home.'data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.gff3';
  	$proteinF = $home.'data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.protein.fa';
	$geneF = $home.'data/IES_data/psexaurelia/sexaurelia_AZ8-4_annotation_v1.gene.fa';
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
    
    my $NF = Bio::SeqIO->new('-file' => $geneF,
			     '-format' => 'Fasta');
    my $PF = Bio::SeqIO->new('-file' => $proteinF,
			     '-format' => 'Fasta');

    print "# reading genes for $species\n";
    while(my $seq = $NF->next_seq()){
	my $id = $seq->primary_id();
	$id =~ s/(P.+)G(\d+)/$1T$2/;
	$nseq{$id} = $seq->seq();
    }
    print "# reading proteins for $species\n";
    while(my $seq = $PF->next_seq()){
	my $id = $seq->primary_id();
	$id =~ s/(P.+)P(\d+)/$1T$2/;
	$pseq{$id} = $seq->seq();
    }

}
#   for each group find how many genes with strange exons and how many total
#loop through group gene allele
#load genes in memory
#load proteins in memory
# my $a = (keys %strangeExons)[0];
# my $b =  (keys %nseq)[0];
# print "$a $b\n";die;
# print "locating 102 genes\n";
open NF, '>'.$home.'data/IES_data/working/genes102N.fa' or die $!;
open PF, '>'.$home.'data/IES_data/working/genes102P.fa' or die $!;
foreach my $gene (keys %strangeExons){
  print $gene,' ',$strangeExons{$gene},"\n";
     if($strangeExons{$gene}>3){
 	if(defined($nseq{$gene})){
 	    print NF '>',$gene, "\n",$nseq{$gene},"\n";
 	}
 	if(defined($pseq{$gene})){
 	    print PF '>',$gene, "\n",$pseq{$gene},"\n";
 	}
     }
}

close NF;
close PF;
# my %countExons;
# my %geneGroups;
# open SLX, $home.'data/IES_data/working/silix.output' or die $!;
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

