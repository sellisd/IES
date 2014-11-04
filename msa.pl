#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
#use Parallel::ForkManager;
#my $cores = 48;
#my $pm = Parallel::ForkManager->new($cores);
#create multiple sequence alignments of silix groups
my $t_coffeeBin = '/panhome/sellis/tools/tcoffee/Version_11.00.8cbe486/bin/t_coffee';
my $dataPath = '/pandata/sellis/';
my $fastaOutPath = '/pandata/sellis/msas/fasta/';
my $logfilePath = '/pandata/sellis/msas/log/';
my $alnfilePath = '/pandata/sellis/msas/aln/';
my $dndfilePath = '/pandata/sellis/msas/dnd/';
my $htmlfilePath = '/pandata/sellis/msas/html/';

#make required directories
mkdir $fastaOutPath unless (-d $fastaOutPath);
mkdir $logfilePath unless (-d $logfilePath);
mkdir $dndfilePath unless (-d $dndfilePath);
mkdir $htmlfilePath unless (-d $htmlfilePath);
mkdir $alnfilePath unless (-d $alnfilePath);
#load groupings in memory
print "loading silix output in memory\n";
my %hash;
my $silixOutput = $dataPath.'working/silix.output';

open IN, $silixOutput or die $!;
while (my $line = <IN>){
    chomp $line;
    (my $group, my $id) = split " ", $line;
    if(defined($hash{$group})){
	push @{$hash{$group}}, $id;
    }else{
	$hash{$group}=[$id];
    }
}
close IN;

#load all protein sequences in memory
print "loading protein sequences in memory\n";
my %proteins;
my @speciesFiles = ('pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.protein.fa','ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.protein.fa','psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.protein.fa');
foreach my $inP (@speciesFiles){
    print "   from $inP\n";
    my $file = Bio::SeqIO->new('-file'=>$dataPath.$inP,
			       '-format' => 'Fasta');
    while(my $seqO = $file->next_seq()){
	my $header = $seqO->primary_id();
	my $sequence = $seqO->seq();
	if(defined($proteins{$header})){
	    die "synonymia?";
	}else{
	    $proteins{$header}  = $sequence;
	}
    }
}

#foreach group make a fasta file with the sequencies
#make list of groups with more than one sequence and align only those
my $groupNo = keys %hash;
print "making fasta files\n";
my @toAlign;
foreach my $group (keys %hash){
    if($#{$hash{$group}}>1){ #make fasta files and align only if thee or more than one sequence in a group
	push @toAlign, $group;
	open OUT, '>'.$fastaOutPath.'/group.'.$group.'.fa' or die $!;
	foreach my $proteinId (@{$hash{$group}}){
	    if(defined($proteins{$proteinId})){
		print OUT ">$proteinId\n$proteins{$proteinId}\n";
	    }else{
		die "$group empty group?";
	    }
	}
	close OUT;
    }
}

#aligning sequences in parallel
print "running t_coffee for each group";
my $dbg = 0;
foreach my $toAlign (@toAlign){
#    my $pid = $pm->start and next;
    #use single core for each group
    my $cmdl = $t_coffeeBin.' '.$fastaOutPath.'group.'.$toAlign.'.fa -multi_core no  &> log.'.$toAlign.'.log';
#make PBS file
    open PBS, '>msa.'.$toAlign.'.pbs' or die $!;
    print PBS '#PBS -q q1hour',"\n";
    print PBS '#PBS -N msa.',$toAlign,"\n";
    print PBS '#PBS -e error.',$toAlign,"\n";
    print PBS '#PBS -o output.',$toAlign,"\n";
    print PBS '#PBS -l nodes=1:dodeca',"\n";
    print PBS "$cmdl","\n";
    print PBS "mv log.*.log $logfilePath\n";
    print PBS "mv group.*.aln $alnfilePath\n";
    print PBS "mv group.*.dnd $dndfilePath\n";
    print PBS "mv group.*.html $htmlfilePath\n";
    print PBS 'echo telos',"\n";
    close PBS;
#    $pm->finish;
    system "qsub msa.".$toAlign.".pbs";
    $dbg++;
    die if $dbg>10;
}
#$pm->wait_all_children;

