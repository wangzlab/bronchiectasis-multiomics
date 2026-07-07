open (IN, $ARGV[0]); ## selected metat and metab module pairs
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$select{$a[0]}{$a[1]}=1;
}
open (IN, "metat_sampleID_match.txt"); ## sample ID match for metat and metab profile for China-GZ cohort
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$match{$a[0]}=$a[1];
}
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

my %metat=();
my %transmetat=();
open (IN, "$sp\_module.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	for my $i (1..$#a) {
		$metat{$a[0]}{$match{$headers[$i]}}=$a[$i];
		$transmetat{$match{$headers[$i]}}{$a[0]}=$a[$i];
	}
}
for my $key2 (keys %select) {
	for my $key3 (keys %{$select{$key2}}) {
		open (OUT, ">metat_metab_mediation_table/$key2\-$key3.txt");
		print OUT "Sample\tMetaT\tMetaB\tPa_Hi\tAge\tGender\tBMI\n";
		for my $key4 (keys %transmetat) {
			#print $key4."\n";
			next unless (exists $transmetab{$key4});
			print OUT $key4."\t".$transmetat{$key4}{$key2}."\t".$transmetab{$key4}{$key3}."\t".$transpahi{$key4}{"Pa_kraken"}."\t".$age{$key4}."\t".$gender{$key4}."\t".$bmi{$key4}."\n";
		}
	}
}
