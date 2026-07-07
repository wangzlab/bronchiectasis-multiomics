open (IN, "Pa_Hi.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	for my $i (1..$#a) {
		$pahi{$a[0]}{$headers[$i]}=$a[$i];
		$transpahi{$headers[$i]}{$a[0]}=$a[$i];
	}
}
open (IN, "gsva.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	for my $i (1..$#a) {
		$hostt{$a[0]}{$headers[$i]}=$a[$i];
		$transhostt{$headers[$i]}{$a[0]}=$a[$i];
	}
}

open (IN, "metab_profile.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	for my $i (1..$#a) {
		$metab{$a[0]}{$headers[$i]}=$a[$i];
		$transmetab{$headers[$i]}{$a[0]}=$a[$i];
	}
}
open (IN, "metadata.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$gender{$a[0]}=$a[2];
	$bmi{$a[0]}=$a[3];
	$age{$a[0]}=$a[4];
	$bsi{$a[0]}=$a[5];
	$smoking{$a[0]}=$a[6];
}
#open (IN, $ARGV[0]);
#while (<IN>) {
#	chop;
#	$select{$_}=1;
#}
for my $key (keys %metab) {
	for my $key2 (keys %hostt) {
		$tmp=$key."-".$key2;
		#next unless (exists $select{$tmp});
		open (OUT, ">metab_hostt_mediation_table/$key\-$key2.txt");
		print OUT "Sample\tMetaB\tHostT\tPa_Hi\tAge\tGender\tBMI\n";
		for my $key3 (keys %transmetab) {
			#print $key4."\n";
			next unless (exists $transhostt{$key3});
			print OUT $key3."\t".$transmetab{$key3}{$key}."\t".$transhostt{$key3}{$key2}."\t".$transpahi{$key3}{"Pa_kraken"}."\t".$age{$key3}."\t".$gender{$key3}."\t".$bmi{$key3}."\n";
		}
	}
}
