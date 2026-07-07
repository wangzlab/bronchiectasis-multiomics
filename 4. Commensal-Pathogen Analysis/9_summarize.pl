open (IN, "metab_desc.txt");
while (<IN>) {
	chop;
	s/ and $//g;
	@a=split("\t",$_);
	$match{"MetaB_ME".$a[0]}=$a[1];
	$metab_desc{"MetaB_ME".$a[0]}=$a[2];
}
open (IN, "hostt_desc.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$match{$a[0]}=$a[1];
	$hostt_desc{$a[0]}=$a[2];
}
open (IN, "KEGG_modules.tab");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$metat_desc{$a[0]}=$a[1];
}
open (IN, "metat_metab_hostt_mediation_results.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$mediate{$a[0]."\t".$match{$a[1]}."\t".$match{$a[2]}}=$a[4]."\t".$a[3]."\t".$a[5]."\t".$a[6]."\t".$a[7]."\t".$a[8]."\t".$a[9]."\t".$a[10]."\t".$a[11]."\t".$a[12]."\t".$a[13]."\t".$a[14];
}
my %direct=();
my %indirect=();
open (IN, $ARGV[0]);
$dump=<IN>;
while (<IN>) {
	chop;
	#s/^map/P/g;
	@a=split("\t",$_);
	next if ($metat_desc{$a[0]} =~ /transport system/);
	$desc{$a[0]."\t".$match{$a[1]}."\t".$match{$a[2]}}=$metat_desc{$a[0]}."\t".$metab_desc{$a[1]}."\t".$hostt_desc{$a[2]};
	$cor{$a[0]."\t".$match{$a[1]}."\t".$match{$a[2]}}=$a[6]."\t".$a[7]."\t".$a[8]."\t".$a[9];
	#$pred{$a[0]."\t".$match{$a[1]}."\t".$match{$a[2]}}=$a[14]."\t".$a[15]."\t".$a[16];
	if (/Directly_Linked/) {
		$direct{$a[0]."\t".$match{$a[1]}."\t".$match{$a[2]}}{$a[10]}=$a[11].", ".$a[13];
	}
	else {
		$indirect{$a[0]."\t".$match{$a[1]}."\t".$match{$a[2]}}{$a[10]}=$a[11].", ".$a[13];

	}
}
print "MetaT\tMetaB\tHostT\tDesc_MetaT\tDesc_MetaB\tDesc_HostT\tCor_MetaT_MetaB\tPval_MetaT_MetaB\tCor_MetaB_HostT\tPval_MetaB_HostT\tDirect\tIndirect\tMetaT_MetaB_FACME\tMetaT_MetaB_FADE\tMetaT_MetaB_FProp\tMetaT_MetaB_RACME\tMetaT_MetaB_RADE\tMetaT_MetaB_RProp\tMetaB_HostT_FACME\tMetaB_HostT_FADE\tMetaB_HostT_FProp\tMetaB_HostT_RACME\tMetaB_HostT_RADE\tMetaB_HostT_RProp\n";
for my $key (sort keys %desc) {
	@tmp=();
	@tmp2=();
	$link=();
	$link2=();
	if (exists $direct{$key}) {
		for my $key2 (keys %{$direct{$key}}) {
			$tmp=$key2." (".$direct{$key}{$key2}.")";
			push @tmp, $tmp;
		}
		$link=join("; ", @tmp);
	}
	else {
		for my $key2 (keys %{$indirect{$key}}) {
			$tmp2=$key2." (".$indirect{$key}{$key2}.")";
			push @tmp2, $tmp2;
		}
		$link2=join("; ", @tmp2);
	}
	print $key."\t".$desc{$key}."\t".$cor{$key}."\t".$link."\t".$link2."\t".$mediate{$key}."\n";
}