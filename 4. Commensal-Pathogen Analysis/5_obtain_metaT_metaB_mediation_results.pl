print "Comparison\tFACME\tFADE\tFTotal_Effect\tFProp_Mediated\tFProp_Mediated_P\tRACME\tRADE\tRTotal_Effect\tRProp_Mediated\tRProp_Mediated_P\n";
open (IN, $ARGV[0]); ## selected metat and metab module list
while (<IN>) {
	chop;
	s/\t/-/g;
	$id=$_;
	open (IN1, "forward/$id\.txt");
	while (<IN1>) {
		chop;
		s/\<//g;
		@a=split(/\s+/,$_);
		if (/ACME/) {
			$facmep=$a[4];
		}
		if (/ADE/) {
			$fadep=$a[4];
		}
		if (/Total/) {
			$ftotalp=$a[4];
		}
		if (/Prop./) {
			$fprop=$a[2];
			$fpropp=$a[5];
		}
	}
	open (IN1, "reverse/$id\.txt");
	while (<IN1>) {
		chop;
		s/\<//g;
		@a=split(/\s+/,$_);
		if (/ACME/) {
			$racmep=$a[4];
		}
		if (/ADE/) {
			$radep=$a[4];
		}
		if (/Total/) {
			$rtotalp=$a[4];
		}
		if (/Prop./) {
			$rprop=$a[2];
			$rpropp=$a[5];
		}
	}
	print $id."\t".$facmep."\t".$fadep."\t".$ftotalp."\t".$fprop."\t".$fpropp."\t".$racmep."\t".$radep."\t".$rtotalp."\t".$rprop."\t".$rpropp."\n";
}
