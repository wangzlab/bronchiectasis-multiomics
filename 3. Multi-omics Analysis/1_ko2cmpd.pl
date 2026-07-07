open (IN, "metacyc_reactions.txt"); ## enzymatic reaction information obtained from MetaCyc
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next if ($a[1] eq "");
	if (/LEFT-TO-RIGHT/) {
		@b=split(";",$a[4]);
		@c=split(";",$a[3]);
		for my $val (@b) {
			$product{$a[1]}{$val}=1;
		}
		for my $val (@c) {
			$substrate{$a[1]}{$val}=1;
		}
	}
	elsif (/RIGHT-TO-LEFT/) {
		@b=split(";",$a[3]);
		@c=split(";",$a[4]);
		for my $val (@b) {
			$product{$a[1]}{$val}=1;
		}
		for my $val (@c) {
			$substrate{$a[1]}{$val}=1;
		}
	}
	elsif (/REVERSIBLE/) {
		$tmp=$a[3].";".$a[4];
		@b=split(";",$tmp);
		for my $val (@b) {
			$product{$a[1]}{$val}=1;
			$substrate{$a[1]}{$val}=1;
		}
	}
}
open (IN, "compounds.tab");  ## compound information obtained from MetaCyc
while (<IN>) {
	chop;
	if (/UNIQUE-ID - (\S+).+HMDB \"HMDB(\d+)\"/) {
		$hmdb{$1}="HMDB00".$2;
	}
}

open (IN, "ko2ec.txt");  ## ko to ec matching information
while (<IN>) {
	chop;
	s/:/-/g;
	@a=split("\t",$_);
	@b=split(" ", $a[1]);
	for my $val (@b) {
		if ($val !~ /EC/) {
			$val = "EC-".$val;
		}
		#print $val."\n";
		$koec{$a[0]}{$val}=1;
	}
}
for my $key (keys %koec) {
	my %koproduct=();
	my %kosubstrate=();
	my %koall=();
	for my $key2 (keys %{$koec{$key}}) {
		next unless (exists $product{$key2} or exists $substrate{$key2});
		for my $key3 (keys %{$product{$key2}}) {
			#if (exists $hmdb{$key3}) {
				#$koproduct{$hmdb{$key3}}=1;
				$koproduct{$key3}=1;
				$koall{$key3}=1;
			#}
		}
		for my $key3 (keys %{$substrate{$key2}}) {
			#if (exists $hmdb{$key3}) {
				#$koproduct{$hmdb{$key3}}=1;
				$kosubstrate{$key3}=1;
				$koall{$key3}=1;
			#}
		}
	}
	#$tmp=join(";", sort keys %koproduct);
	#if ($tmp ne '') {
	#	print $key."\t".$tmp."\n";
	#}
	for my $key4 (keys %koall) {
		if (exists $koproduct{$key4} and exists $kosubstrate{$key4}) {
			print $key."\t".$key4."\t"."Product_And_Substrate\n";
		}
		if (exists $koproduct{$key4} and ! exists $kosubstrate{$key4}) {
			print $key."\t".$key4."\t"."Product\n";		
		}
		if (exists $kosubstrate{$key4} and ! exists $koproduct{$key4}) {
			print $key."\t".$key4."\t"."Substrate\n";
		}
	}
}	
