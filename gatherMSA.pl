#!/usr/bin/perl
use warnings;
use strict;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use lib'.';
use functions;

# gather the data from 3 rounds of MSA on the cluster, if required prepare for local runs
my $path = '/home/dsellis/data/IES/analysis/msas/';
my @runs = qw/1 2/;

foreach my $runNo (@runs){ 
    my $runF = catfile($path,'run'.$runNo.'.dat');
    my $runTF = catfile($path, 'run'.$runNo.'T.dat');
    
    my %map; # map of gene families in runs

    open R, $runF or die $!;
    while(my $line = <R>){
	chomp $line;
	(my $run, my $geneFamily) = split "\t", $line;
	if(defined($map{$run})){
	    push @{$map{$run}}, $geneFamily;
	}else{
	    $map{$run} = [$geneFamily];
	}
    }
    close R;
    
    my $failed = gather({'map'   => \%map,
			 'runTF' => $runTF,
			 'path'  => '/home/dsellis/data/IES/analysis/msas/run'.$runNo,
			 'cpto'  => '/home/dsellis/data/IES/analysis/msas/all/'
			});
    
    print 'round '.$runNo.' had: ', $#{$failed}+1, ' failed runs', "\n";
}

die;
my %failedRunsH;
foreach my $run (@$failed1){
    $failedRunsH{$run} = 1;
}

my $run2F = '/home/dsellis/data/IES/analysis/msas/run2.dat';
my $run2TF = '/home/dsellis/data/IES/analysis/msas/run2T.dat';

my %map2; # map of gene families in runs

open R2, $run2F or die $!;
while(my $line = <R2>){
    chomp $line;
    (my $run2, my $geneFamily) = split "\t", $line;
    if(defined($failedRunsH{$geneFamily})){
	if($failedRunsH{$geneFamily} == 1){
	    # known failure
	}else{
	    die;
	}
    }else{
	die;
    }
    if(defined($map2{$run2})){
	push @{$map2{$run2}}, $geneFamily;
    }else{
	$map2{$run2} = [$geneFamily];
    }
}
close R2;

my $failed2 = gather({'map'   => \%map2,
		  'runTF' => $run2TF,
		  'path'  => '/home/dsellis/data/IES/analysis/msas/run2',
		  'cpto'  => '/home/dsellis/data/IES/analysis/msas/all/'
		 });


print 'round 2 had: ', $#{$failed2}+1, ' failed runs', "\n";


sub gather{
    my $argref = shift @_;
    my $path = $argref->{'path'};
    my $mapref = $argref->{'map'};
    my $cpto = $argref->{'cpto'};
    my @failedRuns;
    make_path($cpto) unless -d $cpto;
    my $alnP = catfile($cpto, 'aln');
    my $htmlP = catfile($cpto, 'html');
    my $dndP = catfile($cpto, 'dnd');
    my $logP = catfile($cpto, 'log');
    make_path($alnP)  unless -d $alnP;
    make_path($htmlP) unless -d $htmlP;
    make_path($dndP)  unless -d $dndP;
    make_path($logP)  unless -d $logP;
    opendir DH, $path or die $!;
    my @dirs = grep {/^cluster.*local$/} readdir(DH);
    open IN, $argref->{'runTF'} or die $!;
    while(my $line = <IN>){
	chomp $line,
	(my $run, my $success) = split "\t", $line;
	if($success == 1){
	    my $pbsF    = catfile($path, 'pbs/msa.'.$run.'.pbs');
	     my $errorF = catfile($path, 'error/error.'.$run);
	    my $outputF = catfile($path, 'output/output.'.$run);
	    if(1 != success({
		'pbs'    => $pbsF,
		'output' => $outputF,
		'error'  => $errorF})){
		die;
	    }
	    # copy to corresponding directories each gene family
	    foreach my $geneFamily(@{$mapref->{$run}}){
		my $found = 0;
		my %files;
		foreach my $node (@dirs){
		    my $alnF = catfile($path, $node, 'cluster.'.$geneFamily.'.aln');
		    my $htmlF = catfile($path, $node, 'cluster.'.$geneFamily.'.html');
		    my $dndF = catfile($path, $node, 'cluster.'.$geneFamily.'.dnd');
		    my $logF = catfile($path, $node, 'cluster.'.$geneFamily.'.fa.log');
		    if(-e $alnF and
		       -e $htmlF and
		       -e $dndF and
		       -e $logF){
			$found++;
			%files = (
			    'alnF'  => $alnF,
			    'htmlF' => $htmlF,
			    'dndF'  => $dndF,
			    'logF'  => $logF,
			);
		    }
		}
		if($found == 1){
		    #copy to new directory
		    my @cmdl = ('cp \''.$files{'alnF'}.'\' '.$alnP,
				'cp \''.$files{'htmlF'}.'\' '.$htmlP,
				'cp \''.$files{'dndF'}.'\' '.$dndP,
				'cp \''.$files{'logF'}.'\' '.$logP);
		    foreach my $cmdl (@cmdl){
#			print $cmdl,"\n";
#			system $cmdl;
		    }
		}elsif($found == 0){
		    print $run," is missing\n";
#		    die;
		}elsif($found >1){
		    print $run, "in multiple nodes\n";
		}else{
		    die;
		}
	    }
	}elsif($success == 0){
	    # build list of gene families for next round
	    push @failedRuns, @{$mapref->{$run}};
	}else{
	    die;
	}	 
    }
    close IN;
    return \@failedRuns;
}

# load all families and search if they are present at least once in one of the run directories
# they need to have 1 pbs file that finished sucessfully

# function sucessful run? pbs file expected output expected error

# build table with results:
# geneFamily inPbs1 run1sucess inPbs2 run2sucess inPbs3 run2sucess
