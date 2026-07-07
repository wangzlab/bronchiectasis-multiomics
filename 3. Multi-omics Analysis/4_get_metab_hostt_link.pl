open (IN, "metabolite_description.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$name{$a[0]}=$a[1];
}
open (IN, "metabo2stitch.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$match{$a[1]}{$a[0]}=1;
}
open (IN, "bronchiectasis_stitch.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next if ($a[2] eq '');
	#next if (/catalysis|reaction/);
	for my $key (keys %{$match{$a[0]}}) {
		$interact{$key}{$a[2]}=$a[4];
		$score{$key}{$a[2]}=$a[9];
	}
}
open (IN, "metab_module_assign.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$metab{"MetaB_ME".$a[1]}{$a[0]}=1;
}
open (IN, "host_pathway_description.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$match2{$a[0]}=$a[1];
}
open (IN, "genego.gmt");
while (<IN>) {
        chop;
        @a=split("\t",$_);
	for my $i (2..$#a) {
	        $trans{$match2{$a[0]}}{$a[$i]}=1;
	}
}
open (IN, $ARGV[0]); ## all_associations.txt from HALLA association between metab modules and host pathways
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next unless ($a[3]<0.05);
	for my $key (keys %{$metab{$a[0]}}) {
		for my $key2 (keys %{$trans{$a[1]}}) {
			if (exists $interact{$key}{$key2}) {
				print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$key."\t".$name{$key}."\t".$key2."\t".$interact{$key}{$key2}."\t".$score{$key}{$key2}."\n";
			}
		}
	}
}

