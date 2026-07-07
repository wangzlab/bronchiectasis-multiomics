open (IN, "metagenome_rel.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	s/ /_/g;
	@a=split("\t",$_);
	$total=0;
	$sum=0;
	for my $i (1..$#a) {
		$metag{$a[0]}{$headers[$i]}=$a[$i];
		$transmetag{$headers[$i]}{$a[0]}=$a[$i];
		$total++;
		$sum+=$a[$i];
	}
	$aver_metag{$a[0]}=$sum/$total;
}
open (IN, "metatranscriptome_rel.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$total=0;
	$sum=0;
	for my $i (1..$#a) {
		$metat{$a[0]}{$headers[$i]}=$a[$i];
		$transmetat{$headers[$i]}{$a[0]}=$a[$i];
		$total++;
		$sum+=$a[$i];
	}
	$aver_metat{$a[0]}=$sum/$total;
}
open (IN, "metadata.txt");
$dump=<IN>;
while (<IN>) {
	chop;
	next unless (/China/);
	next if (/Exacerbations/);
	@a=split("\t",$_);
	if (/Health/) {
		$group{"Healthy"}{$a[0]}=1;
	}
	else {
		$group{$a[1]}{$a[0]}=1;
		#$group{"Stable"}{$a[0]}=1;
	}
}
open (IN, "common_taxa");
while (<IN>) {
	chop;
	$id=$_;
	#print $id."\t".$aver_metag{$id}."\t".$aver_metat{$id}."\n";
	if ($aver_metag{$id}>0.0005 and $aver_metat{$id}>0.0005) {
		$species{$id}=1;
	}
}
for my $key (sort keys %species) {
	$total=0;
	$pos=0;
	$sum=0;
	$sum_simp=0;
	$total2=0;
	$pos2=0;
	$sum2=0;
	$sum_simp2=0;
	$total3=0;
	$pos3=0;
	$sum3=0;
	$sum_simp3=0;
	$total4=0;
	$pos4=0;
	$sum4=0;
	$sum_simp4=0;
	for my $key2 (sort keys %{$group{"Healthy"}}) {
		$total++;
		$sum_simp+=$metag{$key}{$key2};
		if ($metag{$key}{$key2}>0) {
			$pos++;
			$sum+=abs($metag{$key}{$key2}*log($metag{$key}{$key2}));	
		}
	} 
	$aver=$pos/$total * $sum/$total;
	$aver_simp=$sum_simp/$total;
	for my $key2 (sort keys %{$group{"Commensal"}}) {
		$total2++;
		$sum_simp2+=$metag{$key}{$key2};
		if ($metag{$key}{$key2}>0) {
			$pos2++;
			$sum2+=abs($metag{$key}{$key2}*log($metag{$key}{$key2}));	
		}
	}
	$aver2=$pos2/$total2 * $sum2/$total2;
	$aver_simp2=$sum_simp2/$total2;
	$fc=log(($aver+0.000001)/($aver2+0.000001))/log(2);
	$fc_metag{$key}=$fc;
	$fc_metag_simp=log(($aver_simp+0.000001)/($aver_simp2+0.000001))/log(2);
	$fc_metag_simp{$key}=$fc_metag_simp;
	for my $key2 (sort keys %{$group{"Healthy"}}) {
		$total3++;
		$sum_simp3+=$metat{$key}{$key2};
		if ($metat{$key}{$key2}>0) {
			$pos3++;
			$sum3+=abs($metat{$key}{$key2}*log($metat{$key}{$key2}));	
		}
	} 
	$aver3=$pos3/$total3 * $sum3/$total3;
	$aver_simp3=$sum_simp3/$total3;
	for my $key2 (sort keys %{$group{"Commensal"}}) {
		$total4++;
		$sum_simp4+=$metat{$key}{$key2};
		if ($metat{$key}{$key2}>0) {
			$pos4++;
			$sum4+=abs($metat{$key}{$key2}*log($metat{$key}{$key2}));	
		}
	}
	$aver4=$pos4/$total4 * $sum4/$total4;
	$aver_simp4=$sum_simp4/$total4;
	$fc2=log(($aver3+0.000001)/($aver4+0.000001))/log(2);
	$fc_metat{$key}=$fc2;
	$fc_metat_simp=log(($aver_simp3+0.000001)/($aver_simp4+0.000001))/log(2);
	$fc_metat_simp{$key}=$fc_metat_simp;
}
print "Taxa\tAverage_MetaG\tAverage_MetaT\tAM_MetaG\tGM_MetaG\tAM_MetaT\tGM_MetaT\n";
for my $key (sort keys %fc_metag) {
	print $key."\t".$aver_metag{$key}."\t".$aver_metat{$key}."\t".$fc_metag_simp{$key}."\t".$fc_metag{$key}."\t".$fc_metat_simp{$key}."\t".$fc_metat{$key}."\n";
}

