sub buildPaths{
    # return standard file names and paths
    my $binomial = shift @_;
    my $pabAnot = shift @_;
    $binomial =~ /^([A-Z])[a-z]+\s+([a-z]+)$/ or die;
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

sub prefix{
    my $species = shift @_;
    my %prefixes = (
	'ppr' => 'PPRIM.AZ9-3.1.',
	'pbi' => 'PBIA.V1_4.1.',
	'pte' => 'PTET.51.1.',
	'ppe' => 'PPENT.87.1.',
	'pse' => 'PSEX.AZ8_4.1.',
	'poc' => 'POCTA.138.1.',
	'ptr' => 'PTRED.209.2.',
	'pca' => 'PCAU.43c3d.1.'
#	'tth' => ''	   
    );
    if (defined($prefixes{$species})){
	return $prefixes{$species};
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

# sub mergeOverlap{
# 	# if elements in @A partially overlap return which should be merged
# 	my $asref = shift @_;
# 	my $aeref = shift @_;
#     my $n = $#{$asref} +1;
#     for (my $i = 0; $i <$n; $i++){
# 		for (my $j = 0; $j <$n; $j++){
# 			my $overlap = 0;
# 			if(@$asref[$i]>=@$asref[$j] and @$asref[$i] <= @$aeref[$j]){
# 				$overlap = 1;
# 			}
# 			if(@$aeref[$i]>=@$asref[$j] and @$aeref[$i] <= @$aeref[$j]){
# 				$overlap = 1;
# 			}
# 		}
# 	}
# 	# group pairs

# }
# 1;

# 0 --- 
# 1   ---
# 2    ----

# starti within j
# or
# endi within j

# 0.1 =>{0=>1,
# 	   1=>1}
# 0.2 => {0=>1,
#         2=>1}
# 1.2 => {1=>1,
#         2=>1}

# 0.1=> {0=>}
1;
