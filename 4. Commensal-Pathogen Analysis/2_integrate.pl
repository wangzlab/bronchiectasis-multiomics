open (IN, $ARGV[0]); ## metat_metab link for each selected species
$dump=<IN>;
while (<IN>) {
	chop;
	next if (/Same_MetaG_MetaB_Module/);
	@a=split("\t",$_);
	next unless ($a[11]<0.1);
	$mode{$a[1]."\t".$a[0]}{$a[4]."\t".$a[5]}=$a[8];
	$correlation{$a[1]."\t".$a[0]}=$a[10]."\t".$a[11];
	$desc{$a[0]}=$a[2];
	$desc{$a[1]}=$a[3];
	$link{$a[1]."\t".$a[0]}{$a[4]."\t".$a[5]}=1;
}

open (IN, "host_pathway_description.txt"); 
while (<IN>) {
	chop;
	#s/results_/HostT_ME/g;
	@a=split("\t",$_);
	$desc{$a[1]}=$a[0];
}
open (IN, "metab_hostt_link.txt"); ## metab host link select
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next unless ($a[3]<0.1);
	$correlation2{$a[0]."\t".$a[1]}=$a[2]."\t".$a[3];
	$link2{$a[0]."\t".$a[1]}{$a[4]."\t".$a[6]."\t".$a[7]}=1;
}
open (IN, "metabolite_description.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$desc{$a[0]}=$a[1];
}
print "MetaG\tMetaB\tHostT\tDesc_MetaG\tDesc_MetaB\tDesc_HostT\tCor_MetaG_MetaB\tPval_MetaG_MetaB\tCor_MetaB_HostT\tPval_MetaB_HostT\tLinks\tMetabolite_Description\tWays_of_Link1\tWays_of_Link2\n";
for my $key (keys %link) {
	($metag,$metab)=($key =~ /(\S+)\t(MetaB_\S+)/);
	for my $key2 (keys %link2) {
		($metab2,$trans)=($key2 =~ /(MetaB_\S+)\t(\S+)/);
		next unless ($metab2 eq $metab);
		#print $metab2."\t".$metab."\n";
		for my $key3 (keys %{$link{$key}}) {
			($metagenome,$metabolite)=($key3 =~ /(\S+)\t(\S+)/);
			for my $key4 (keys %{$link2{$key2}}) {
				($metabolite2,$hostgene,$int)=($key4 =~ /(\S+)\t(\S+)\t(\S+)/);
				next if ($int =~ /catalysis|reaction/);
				#print $metabolite2."\t".$hostgene."\t".$int."\n";
				if ($metabolite eq $metabolite2) {
					print $metag."\t".$metab."\t".$trans."\t".$desc{$metag}."\t".$desc{$metab}."\t".$desc{$trans}."\t".$correlation{$key}."\t".$correlation2{$key2}."\t".$metagenome."_".$metabolite."_".$hostgene."\t".$desc{$metabolite}."\t".$mode{$key}{$key3}."\t".$int."\n";
				}
			}
		}
	}
}
