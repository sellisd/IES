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
