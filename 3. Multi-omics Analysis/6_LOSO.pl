open (IN, "bacterial_species_code.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$sp{$a[0]}{$a[1]}=1;
	$sp2{$a[1]}=$a[0];
}
open (IN, "ko_noeukaryotes.txt");
while (<IN>) {
        chop;
        @a=split("\t",$_);
	($species)=($a[1] =~ /(\S+)\:/);
	next unless (exists $sp2{$species});
	$ko{$sp2{$species}}{$a[0]}=1;
}
$sp=$ARGV[0];
open (IN, "all_gene_depth.txt");
open (OUT, ">all_gene_depth_$sp.txt");
$header=<IN>;
print OUT $header;
#chop $header;
#@headers=split("\t",$header);
while (<IN>) {
        chop;
        @a=split("\t",$_);
	next if (exists $ko{$sp}{$a[0]});
	print OUT $_."\n";
}

