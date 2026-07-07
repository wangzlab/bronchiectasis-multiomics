open (IN, "cmpd2metabo.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$hash{$a[1]}{$a[0]}=1;
}
open (IN, "ko2cmpd.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next unless (exists $hash{$a[1]});
	for my $key (keys %{$hash{$a[1]}}) {
		#print $a[1]."\t".$key."\n";
		if ($a[2] eq 'Product') {
			$product{$a[0]}{$key}=1;
			$all{$a[0]}{$key}=1;
		}
		if ($a[2] eq 'Substrate') {
			$substrate{$a[0]}{$key}=1;
			$all{$a[0]}{$key}=1;
		}
		if ($a[2] eq 'Product_And_Substrate') {
			$product{$a[0]}{$key}=1;
			$substrate{$a[0]}{$key}=1;
			$all{$a[0]}{$key}=1;
		}
	}
}
for my $key (sort keys %all) {
	for my $key2 (sort keys %{$all{$key}}){
		if (exists $product{$key}{$key2} and exists $substrate{$key}{$key2}) {
			print $key."\t".$key2."\tProduct_And_Substrate\n";
		}
		if (exists $product{$key}{$key2} and ! exists $substrate{$key}{$key2}) {
			print $key."\t".$key2."\tProduct\n";
		}

		if (! exists $product{$key}{$key2} and exists $substrate{$key}{$key2}) {
			print $key."\t".$key2."\tSubstrate\n";
		}
	}
}
