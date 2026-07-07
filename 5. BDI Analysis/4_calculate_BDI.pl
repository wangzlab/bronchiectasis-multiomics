open (IN, "combined_metagenome_three_controls.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	s/\S+s__/s__/g;
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
open (IN, "fold_change_metag_metat.txt");
while (<IN>) {
	chop;
	s/\S+s__/s__/g;
	@a=split("\t",$_);
	$fc{$a[0]}=$a[1];
}
#for (my $i=0;$i<=3;$i+=0.05) {  # threshold for selecting health-related FC cutoff
#	for (my $j=0;$j<=5;$j+=0.05) { # threshold for selecting disease-related FC cutoff
	for $i (0.85) {
	for my $j (1.1) {
		%disease=();
		%health=();
		open (IN, "pathogen_list");
		while (<IN>) {
			chop;
			s/\S+s__/s__/g;
			$disease{$_}=1;
		}
		$m=0;
		$n=0;
		for my $key (keys %fc) {
			if ($fc{$key}>=$i) {
				$health{$key}=1;
				$m++;
				#print STDERR $i."\t".$j."\t".$key."\t"."Health\n";
			}
			if ($fc{$key}<=-1*$j) {
				$disease{$key}=1;
				$n++;
				#print STDERR $i."\t".$j."\t".$key."\t"."Disease\n";
			}
		}
		next if ($m==0 or $n==0);
		$h=0;
		$hh=0;
		$d=0;
		$dd=0;
		#open (OUT, ">CI_diff_cutoffs_AM_old2/CI_${i}_${j}.txt");
		for my $key (sort keys %transmetag) {
			$totalY=0;
			$posY=0;
			$sumY=0;
			for my $key2 (keys %health) {
				$totalY++;
				if ($metag{$key2}{$key}>0) {
					$posY++;
					$tmp=$metag{$key2}{$key};
					#$sumY+=abs($tmp*log($tmp));
					#$sumY+=log($tmp+0.000001);
					$sumY+=log($tmp);
				}
			}
			#$phiY=$posY/$totalY*$sumY;
			#$phiY=$sumY/$totalY;
			$totalN=0;
			$posN=0;
			$sumN=0;
			for my $key2 (keys %disease) {
				$totalN++;
				if ($metag{$key2}{$key}>0) {
					$posN++;
					$tmp=$metag{$key2}{$key};
					#$sumN+=abs($tmp*log($tmp));
					#$sumN+=log($tmp+0.000001);
					$sumN+=log($tmp);
				}
			}
			next if ($posN==0 or $posY==0);
			$phiY=$sumY/$posY;
			#$phiN=$sumN/$totalN;
			$phiN=$sumN/$posN;
			#$phiN=$posN/$totalN*$sumN;
			#$diff=log(($phiY+0.000001)/($phiN+0.000001))/log(10);
			$diff=$phiN-$phiY;
			print $key."\t".$phiY."\t".$phiN."\t".$diff."\n";
		}
	}
}
