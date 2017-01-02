#!/usr/bin/perl
use File::HomeDir;
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path);
use lib'.';
use functions;
use warnings;
use strict;

# Evolution of IES in *Paramecium*

# load options and variables
my $opt = loadUserOptions;
my $basePath = $$opt{'basePath'};
my $notationF =  catfile($basePath, 'analysis/notation.csv');
my $nr = getNotation($notationF);
my %notation;

# Variables
my $scafCutoff = $$opt{'scafCutoff'};
my $tablesP    = catdir($basePath, 'analysis', 'tables');
my $bedP       = catdir($basePath, 'analysis', 'bed');
my $filtscafP  = catdir($basePath, 'analysis', 'filtscaf');
my $iesdbP     = catdir($basePath, 'analysis', 'iesdb');
my $protdbP = catdir($basePath, 'analysis', 'protdb');
make_path($protdbP)   unless -d $protdbP;
make_path($bedP)      unless -d $bedP;
make_path($tablesP)   unless -d $tablesP;
make_path($filtscafP) unless -d $filtscafP;
make_path($tablesP)   unless -d $tablesP;

# calculate scaffold lengths
# --------------------------
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}}; #de-reference for less typing
    my $genome = catfile($pab{'datapath'}, $pab{'MacF'});
    my $scafF = catfile($basePath, 'analysis', 'filtscaf', $pab{'abr'}.'.scaf');
    my $cmdl = './scaffoldStats.pl '.$genome.' > '.$scafF;
    run($cmdl, 1);
}

# filter scaffolds
# ----------------

foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $scafF   = catfile($filtscafP, $pab{'abr'}.'.scaf');
    my $protein = catfile($pab{'datapath'}, $pab{'protF'});
    my $gff     = catfile($pab{'datapath'}, $pab{'geneGff'});
    my $gene    = catfile($pab{'datapath'}, $pab{'geneF'});
    my $ies     = catfile($pab{'datapath'}, $pab{'iesGff'});
    my $outdir  = $filtscafP;
    my $cutoff  = $scafCutoff;
    my $cmdl = './filterScaffolds.pl'.
               ' -length  '.$scafF.
               ' -ies     '.$ies.
               ' -protein '.$protein.
               ' -gene    '.$gene.
               ' -gff     '.$gff.
               ' -cutoff  '.$cutoff.
               ' -species '.$pab{'abr'}.
               ' -outdir  '.$outdir;
    run($cmdl, 1);
}


# make BLAST protein database
# ---------------------------

#find /home/dsellis/data/IES/analysis/ -name "*.protein.fa"

my $allprotF = catfile($protdbP, 'allprot.fa');
my $tthGenome = catfile($basePath, 'genomicData', 'thermophila', 'gene', 'tth.protein.fa');
my @proteinF;
for my $sp (keys %$nr){
    my %pab = %{$nr->{$sp}};
    push @proteinF, catfile($filtscafP, $pab{'abr'}.'.protein.fa');
}
push @proteinF, $tthGenome; #add T. thermophila
my $cmdl = 'cat '.join(' ', @proteinF).' > '.$allprotF;
run($cmdl, 1);
$cmdl = 'makeblastdb -in '.$allprotF.' -dbtype prot -out '.catfile($protdbP,'allprot');
run($cmdl, 1);

# split fasta files for faster processing in parallel
# ---------------------------------------------------
my $allprotFa = catfile($basePath, 'analysis', 'protdb', 'allprot.fa');
my $fastaChunks = catfile($basePath, 'tempdat', 'fastachunks');
run('./preblast.pl '.$allprotFa.' '.$fastaChunks, 1);

# parse basic ies information
#----------------------------
for my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $gffF = catfile($filtscafP, $pab{'abr'}.'.ies.gff3');
    my $tabF = catfile($filtscafP, $pab{'abr'}.'.ies.tab');
    my $cmdl = "./iesInfo.pl $gffF > $tabF";
    run($cmdl, 1);
}

# find which IES are floating
for my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $tabF   = catfile($filtscafP, $pab{'abr'}.'.ies.tab');
    my $floatF = catfile($filtscafP, $pab{'abr'}.'.ies.float');
    my $bedF   = catfile($bedP, $pab{'abr'}.'.ies.be'); 
    my $cmdl = "./annotateFloating.pl -tabF $tabF -floatF $floatF -bedF $bedF";
    run($cmdl, 1);
}

# check if we need to merge floating IES
for my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $gffF = catfile($filtscafP, $pab{'abr'}.'.gff');
    my $bedF = catfile($bedP, $pab{'abr'});
    my $cmdl = "./geneInfo.pl $gffF $bedF";
    run($cmdl, 1);
}

# create bed files
for my $sp (sort keys %$nr){
    my %pab     = %{$nr->{$sp}};
    my $exonF   = catfile($bedP, $pab{'abr'}.'.exon.be');
    my $intronF = catfile($bedP, $pab{'abr'}.'.intron.be');
    my $cmdl = "Rscript --vanilla exon2intron.R $exonF $intronF";
    run($cmdl, 1);
}

for my $sp (sort keys %$nr){
    my %pab     = %{$nr->{$sp}};
    my $geneF   = catfile($bedP, $pab{'abr'}.'.gene.be');
    my $interF  = catfile($bedP, $pab{'abr'}.'.inter.be');
    my $scafF   = catfile($basePath, 'analysis', $pab{'abr'}.'.scaf');
    my $cmdl = "Rscript --vanilla gene2intergenic.R $geneF $interF $scafF";
    run($cmdl, 1);
}

# find overlap in bed files
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $iesF    = catfile($bedP, $abr.'.ies.be');
    my $cdsF    = catfile($bedP, $abr.'.cds.be');
    my $intronF = catfile($bedP, $abr.'.intron.be');
    my $interF  = catfile($bedP, $abr.'.inter.be');
    my $iesinF  = catfile($bedP, $abr.'.IESin.be'); 
    my $cmdl = "bedtools intersect -a $iesF".
	" -b $cdsF".
	" -b $intronF".
	" -b $interF".
	" -names cds intron intergenic -wo > $iesinF";
    run($cmdl, 1);
}

# make IES table for iesDB
my $mergeF = catfile($filtscafP, 'ies2merge.dat');
$cmdl = "./reannotateFloating.pl ".$filtscafP."/*.ies.float > $mergeF";
run($cmdl, 1);

foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $iesinF  = catfile($bedP, $abr.'.IESin.be');
    my $iestabF = catfile($filtscafP, $abr.'.ies.tab');
    my $floatF  = catfile($filtscafP, $abr.'.ies.float');
    my $iesdbF  = catfile($iesdbP, $abr.'.iesdb');
    my $cmdl = './iesdbTable.pl -merge '.$mergeF.
	" -iesin  $iesinF".
	" -iestab $iestabF".
    " -float  $floatF".
	" > $iesdbF";
    run($cmdl, 1);
}

# make CDS table for iesDB
my $silixout = catfile($basePath, 'analysis', 'allvsall', 'blastout', 'silix.output');
foreach my $sp (sort keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $cdsdbF = catfile($iesdbP, $abr.'.cdsdb');
    my $cmdl = './cdsdbTable.pl -silixout '.$silixout.
	' -gff '.catfile($pab{'datapath'}, $pab{'geneGff'}).
	' > '.$cdsdbF;
    run($cmdl, 1);
}

# find transcript coordinates for IES in genes
foreach my $sp (keys %$nr){
    my %pab = %{$nr->{$sp}};
    my $abr = $pab{'abr'};
    my $iesinF      = catfile($bedP, $abr.'.IESin.be');
    my $cdsdbF      = catfile($iesdbP, $abr.'.cdsdb');
    my $iesInGenesF = catfile($tablesP, $abr.'.iesInGenes');
    my $cmdl = 'Rscript --vanilla ./iesInGenes.R '.$iesinF.' '.$cdsdbF.' '.$iesInGenesF;
    run($cmdl, 1);
}
