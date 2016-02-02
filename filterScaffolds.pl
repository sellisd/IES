#!/usr/bin/perl
use warnings;
use strict;
use lib'.';
use File::Path qw(make_path);
use functions;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
use File::Spec::Functions qw(catfile);
my $help;
my $cutoff = 10**4;
my $lengthF;
my $iesF;
my $geneF;
my $proteinF;
my $gff;
my $species;
my $outdir;
my $usage = <<HERE;

Filter proteins, genes and IES in scaffolds larger than a size cutoff
usage: filterSCaffolds.pl [OPTIONS]
where OPTIONS can be:
    -length   : table with scaffold lengths (output of scaffoldStats.pl)
    -ies      : path and file name of IES file (GFF3)
    -protein  : path and file name of protein file (fasta)
    -gene     : path and file name of gene file (fasta)
    -gff      : main annotation file
    -cutoff   : cuttoff, default 10^4
    -species  : species three letter abreviation
    -outdir   : directory with output filtered files
    -help|?: this help screen

HERE

die $usage unless (GetOptions(
		       'help|?'     => \$help,
		       'gff=s'      => \$gff,
		       'cutoff=i'   => \$cutoff,
		       'length=s'   => \$lengthF,
		       'protein=s'  => \$proteinF,
		       'gene=s'     => \$geneF,
		       'ies=s'      => \$iesF,
		       'species=s'  => \$species,
		       'outdir=s'   => \$outdir
		   ));

die $usage if $help;
my $prefix = &prefix($species);
die unless $prefix;
make_path($outdir) unless -d $outdir;

# find which scaffolds are larger than the cutoff
print "Read scaffold lengths...";
my %scL;
open L, $lengthF or die "$lengthF $!";
my $header = readline(L);
my $counter;
my $totalC;
while(my $line = <L>){
    chomp $line;
    (my $scaffold, my $length) = (split " ", $line)[0,1];
    if ($length > $cutoff){
	$scL{$scaffold} = $length;
	$counter++;
    }
    $totalC++;
}
close L;
print "done\n";
print "$counter out of $totalC scaffolds over the cutoff\n";

# find which genes are in filtered scaffolds
print "Building scaffold/gene hash...";
my %genesInScaf;
my $gffO = Bio::Tools::GFF->new('-file'        => $gff,
				'-gff_version' => 3);
my $geneCount;
my $geneCountFilt;
while(my $feature = $gffO->next_feature()){
    my $type = $feature->primary_tag;
    next unless($type eq 'gene');
    $geneCount++;
    my $name = $feature->primary_id;
    $name =~ s/($prefix)G// or die;
    my $scaffold = $feature->seq_id;
    if(defined($scL{$scaffold})){
	$genesInScaf{$name} = $feature->seq_id;
	$geneCountFilt++;
    }
   #  print $feature->primary_id, "\n";
   #  print $feature->primary_tag, "\n";
   #  print $feature->source_tag, "\n";
   #  print $feature->seq_id, "\n";
}
print "done\n";
print "$geneCountFilt out of $geneCount genes over the cutoff\n";

#read proteins and filter
my $protIn = Bio::SeqIO->new('-file'   => $proteinF,
			     '-format' => 'fasta');
my $protoutF = catfile($outdir, $species.'.protein.fa');
my $protOut = Bio::SeqIO->new('-file'  => '>'.$protoutF,
			      '-format' => 'fasta');
while(my $seqO = $protIn->next_seq){
    my $protName =  $seqO->primary_id;
    $protName =~ s/($prefix)P// or die;
    if(defined($genesInScaf{$protName})){
	$protOut->write_seq($seqO);
    }
}

#read genes and filter
my $geneIn = Bio::SeqIO->new('-file'    => $geneF,
			     '-format'  => 'fasta');
my $geneoutF = catfile($outdir, $species.'.gene.fa');
my $geneOut = Bio::SeqIO->new('-file'   => '>'.$geneoutF,
			      '-format' => 'fasta');

while(my $seqO = $geneIn->next_seq){
    my $geneName =  $seqO->primary_id;
    $geneName =~ s/($prefix)G// or die;
    if(defined($genesInScaf{$geneName})){
	$geneOut->write_seq($seqO);
    }
}

#read IES and filter
my $iesIn = Bio::Tools::GFF->new('-file'         => $iesF,
				 '-gff_version'  => 3);
my $outiesF = catfile($outdir, $species.'.ies.gff3');
my $iesOut = Bio::Tools::GFF->new('-file'        => '>'.$outiesF,
				  '-gff_version' => 3);
my $iesCount;
my $iesCountFilt;
while(my $feature = $iesIn->next_feature()){
    my $scaffold = $feature->seq_id;
    $iesCount++;
    if(defined($scL{$scaffold})){
	$iesOut->write_feature($feature);
	$iesCountFilt++;
    }
}
print "$iesCountFilt out of $iesCount IES over the cutoff\n";

