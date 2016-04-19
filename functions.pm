package functions;
use strict;
use warnings;
BEGIN{
    require Exporter;
    our $VERSION = 1.01;
    our @ISA = qw(Exporter);
    our @EXPORT = qw(run getNotation isFloating prot2gene printab buildPaths abr2prefix whichInOne gene2species gene2prot initF success);
    our @EXPORT_OK = qw();
}

sub run{
    my $cmdl = shift @_;
    my $dryRun = shift @_;
    print $cmdl,"\n";
    system $cmdl unless($dryRun);
}

sub getNotation{
    my $notationF = shift @_; # file with standard notation
    open N, $notationF or die $!;
    my $header = readline(N);
    my %notation;
    while(my $line = <N>){
	chomp $line;
	(my $abbreviation, my $datapath, my $binomial, my $taxId, my $geneGff, my $cdsF, my $protF, my $geneF, my $MacF, my $iesGff, my $annotation, my $prefix) = split "\t", $line;
	$notation{$binomial} = {
	    'annotation' => $annotation,
	    'prefix'     => $prefix,
	    'abr'        => $abbreviation,
	    'datapath'   => $datapath,
	    'taxId'      => $taxId,
	    'geneGff'    => $geneGff,
	    'cdsF'       => $cdsF,
	    'protF'      => $protF,
	    'geneF'      => $geneF,
	    'MacF'       => $MacF,
	    'iesGff'     => $iesGff
	};
    }
    close N;
    return \%notation;
}

sub success{
    # check if a pbs run was successful
    # by having an error file with zero size or with an expected warning we can ignore
    my $argref = shift @_;
    my $pbsF = $argref->{'pbs'};
    my $outputF = $argref->{'output'};
    my $errorF = $argref->{'error'};
    my $success = 0;
    if(! -e $errorF){   # does not exist
	$success = 0;
    }else{
	if(-z $errorF){     # has zero size
	    $success = 1;
	}else{              # has non-zero size
	    open ERR, $errorF or die $!;
	    my $line = readline(ERR);
	chomp $line;
	    if ($line =~ /^mkdir: cannot create directory(.*)File exists$/){
		$success = 1;
	    }
	    close ERR;
	}
    }
    return $success;
}

sub initF{
    my $notationF =  '/home/dsellis/data/IES/analysis/notation.tab';
    open N, $notationF or die $!;
    my $header = readline(N);
    my %prefixes;
    while(my $line = <N>){
	chomp $line;
	my $ar = split "\t", $line;
	(my $abbreviation, my $datapath, my $binomial, my $taxId, my $geneGff, my $cdsF, my $protF, my $geneF, my $MacF, my $iesGff, my $annotation, my $prefix) = split "\t", $line;
	$prefixes{$abbreviation} = $prefix;
    }
    close N;
    return \%prefixes;
}


# my %prefixes = (
#     'ppr' => 'PPRIM.AZ9-3.1.',
#     'pbi' => 'PBIA.V1_4.1.',
#     'pte' => 'PTET.51.1.',
#     'ppe' => 'PPENT.87.1.',
#     'pse' => 'PSEX.AZ8_4.1.',
#     'poc' => 'POCTA.138.1.',
#     'ptr' => 'PTRED.209.2.',
#     'pca' => 'PCAU.43c3d.1.'
# #	'tth' => ''	   
#     );

sub isFloating{
    # check if an IES is floating
    # returns 0 for not floating, +? or -? if flank sequence is too short and an array of displacements in case of floating
    my $iesS = shift @_;
    my $macUpstreamS = shift @_;
    my $macDownstreamS = shift @_;
    # validate input
    die "input sequence error:\n",
        "  IES:              $iesS\n",
        "  upstream flank:   $macUpstreamS\n",
        "  downstream flanK: $macDownstreamS" unless($iesS =~ /^TA[ACTGN]+TA$/
						    and $macUpstreamS =~ /^[ACTGN]+TA$/
						    and $macDownstreamS =~ /^TA[ACTGN]+$/);
    # transform strings to arrays for easy handling 
    my @iesA = split('', $iesS);
    my @macUpstreamA = split('', $macUpstreamS);
    my @macDownstreamA = split('', $macDownstreamS);
    
    #iesS is an array of IES sequence including both TAs
    #macUpstreamS and macDownstreamS is an array both including the TA junction
    
    #default values
    my @altLoc;
    my $floating = 0;
    my $altLoc;
    #if sequence starts with TA(x)kTA... and downstream Mac boundary is (x)kTA
    #or if sequence ends with ...TA(x)kTA and upstream Mac boundary is TA(x)k
    my $maxSeq = length($iesS);
    # search downstream first
    for(my $i = 2; $i < $maxSeq; $i++){
        if($i > $#macDownstreamA){
            return '+?'; # Did not reach end of downstream comparison
        }
#      print "$i $iesA[$i] $macDownstreamA[$i] ";
	if($iesA[$i] eq $macDownstreamA[$i]){
	    # keep searching in this direction
#	    print "  keep going downstream\n";
	    if($iesA[$i - 1] eq 'T' and $iesA[$i] eq 'A'){
		# is floating
		$floating = 1;
		$altLoc = $i - 1; # location of last matching TA
		push @altLoc, $altLoc;
	    }       
	}else{
	    last;
	}
    }
    $floating = 0;
    # search upstream then
    for(my $i = 2 ; $i < $maxSeq; $i++){
	if($i > $#macUpstreamA){
	    return '-?'; # Did not reach end of upstream comparison
	}
	if($iesA[$#iesA - $i] eq $macUpstreamA[$#macUpstreamA - $i]){ # compare from end
	    # keep searching in this direction
	    if($iesA[$#iesA - $i] eq 'T' and $iesA[$#iesA  - $i + 1] eq 'A'){
		# is floating
		$floating = 1;
		$altLoc = -($i - 1); # location of last matching TA
		push @altLoc, $altLoc;
	    }
	}else{
	    last;
	}
    }
    if(@altLoc){
	return \@altLoc;
    }else{
	return 0;
    }
# for i++
# compare ies begin/end to mac up/down
# if identical continue
# if not was the previous one a T-T match or an A-A match and currently the mac has an A or a T?
# floating
# if end of ies and fully identical this is a potential repeated floating
    
}

sub prot2gene{
    my $protId = shift @_;
    my $prefixesR = shift @_;
    my $found = 0;
    my $geneId;
    foreach my $abr (keys %$prefixesR){
	my $prefix = $prefixesR->{$abr};
	if($protId =~ /^$prefix[P](\d+)$/){
	    $found++;
	    $geneId = $prefix.'G'.$1;
	}
    }
    if($found == 1){
	return $geneId;
    }else{
	return;
    }  
}

sub gene2prot{
    my $geneId = shift @_;
    my $prefixesR = shift @_;
    my $found = 0;
    my $protId;
    foreach my $abr (keys %$prefixesR){
	my $prefix = $prefixesR->{$abr};
	if($geneId =~ /^$prefix[G](\d+)$/){
	    $found++;
	    $protId = $prefix.'P'.$1;
	}
    }
    if($found == 1){
	return $protId;
    }else{
	return;
    }
}

sub printab{
    print(join("\t", @_),"\n");
}

sub buildPaths{
    # return standard file names and paths
    my $binomial = shift @_;
    my $pabAnot = shift @_;
    $binomial =~ /^([A-Z])[a-z]+\s+([a-z]+)$/ or die $binomial;
    my $g = lc($1);
    my $spEpithet = $2;
    my $pab = $g.substr($spEpithet, 0, 2);
    my $scaffoldF = 'data/IES/analysis/'.$pab.'.scaf';
    my $proteinF  = 'data/IES/'.$spEpithet.'/gene/'.$pabAnot.'.protein.fa';
    my $geneF     = 'data/IES/'.$spEpithet.'/gene/'.$pabAnot.'.gene.fa';
    my $anotF     = 'data/IES/'.$spEpithet.'/gene/'.$pabAnot.'.gff3';
    my $iesF      = 'data/IES/'.$spEpithet.'/IES/'.$g.$spEpithet.'_internal_eliminated_sequence.gff3';
    return {
	'scaffoldF' => $scaffoldF,
	'proteinF'  => $proteinF,
	'geneF'     => $geneF,
	'anotF'     => $anotF,
	'iesF'      => $iesF,
	'pab'       => $pab
    }
}

sub abr2prefix{
    my $abr = shift @_;
    my $prefixesR = shift @_;
    die "missing prefixes" unless defined($prefixesR);
    if (defined($prefixesR->{$abr})){
	return $prefixesR->{$abr};
    }else{
	return;
    }
}

sub whichInOne{
# find which elements in @A are completely within ONE element of @B
    my $asref = shift @_;
    my $aeref = shift @_;
    my $bsref = shift @_;
    my $beref = shift @_;
    my %inB;
    my $na = $#{$asref} + 1;
    my $nb = $#{$bsref} + 1;
    for (my $i = 0; $i <$na; $i++){
		for (my $j = 0; $j <$nb; $j++){
	    	if(@$asref[$i] >= @$bsref[$j]){
				if(@$aeref[$i] <= @$beref[$j]){
		    		# completely inside
		    		if(defined($inB{$i})){
						$inB{i} = 0;
		    		}else{
						$inB{$i} = 1;
			    	}
				}
	    	}
		}
    }
    my %inOne;
    foreach my $entry (keys %inB){
		if ($inB{$entry} == 1){
	    	$inOne{$entry} = 1;
		}
	}
    return \%inOne;
}

sub gene2species{
#get species name from gene name
    my $string = shift @_;
    my $speciesName;
    my $abr = substr($string,0,4);
    if($abr eq 'PCAU'){
	$speciesName = 'Paramecium_caudatum';
    }elsif($abr eq 'PSEX'){
	$speciesName = 'Paramecium_sexaurelia';
    }elsif($abr eq 'PTET'){
	$speciesName = 'Paramecium_tetraurelia';
    }elsif($abr eq 'PBIA'){
	$speciesName = 'Paramecium_biaurelia';
    }elsif($abr eq 'TTHE'){
	$speciesName = 'Tetrahymena_thermophila';
    }else{
	die "unknown name $abr";
    }
}

1;
