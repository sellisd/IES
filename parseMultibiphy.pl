#!/usr/bin/perl
use warnings;
use strict;

my $burnIn = 100;
my $mbfPath = '/home/dsellis/data/IES_data/msas/asr/results/';

opendir DH, $mbfPath or die $!;
my @files = grep {/^ies\.cluster\.\d+\.phy\.treelist$/} readdir(DH);
close DH;
foreach my $fileName (@files){
#    $fileName = 'ies.cluster.10.phy.treelist';
    open IN, $mbfPath.$fileName or die $!;
    my $lineCounter = 0;
    my $totalCycles = 1;
    my @counts;
    while(my $line = <IN>){
	chomp $line;
	my @ar = split " ", $line;
	my $iesNo; #number of IES in cluster
	if($lineCounter == 0){ #header
	    my @nodes = @ar;
	}else{
	    if($lineCounter > $burnIn){
		#in each column split by character
		my $colCounter = 0;
		foreach my $entry (@ar){
		    my @presenceAbsence = split "", $entry;
		    $iesNo = $#presenceAbsence + 1;
		    for(my $i = 0; $i <= $#presenceAbsence; $i++){
			if($presenceAbsence[$i] == 1){
			    if(defined($counts[$colCounter+$i])){
				$counts[2*$colCounter+$i]++;
			    }else{
				$counts[2*$colCounter+$i] = 1;
			    }
			}
		    }
		    $colCounter++;
		}
		$totalCycles++;
	    }
#array size is number of IES
	    #{node => [ies1Counter,ies2Counter,...]}
	}
	$lineCounter++;
    }
    close IN;
    $fileName =~ /^ies\.cluster\.(\d+)\.phy\.treelist$/;
    my $cluster = $1;
    print $cluster,"\t";
    foreach my $sum (@counts){
	if(defined($sum)){
	    print $sum/($totalCycles-1),"\t";
	}else{
	    print "0\t";
	}
    }
    print "\n";
    #print $totalCycles-1,"\n"; 
#    die;
}

# for branch Si to Sj
# calculate difference from all Si to Sj (including paralogs)

# for all pair of paralogs compare Si to Sj (how similar they are)
