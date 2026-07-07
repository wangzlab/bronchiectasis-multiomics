$sp=$ARGV[0];
open (IN, "$sp\.txt"); ## KO file for each selected commensal species
$dump=<IN>;
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$ko{$a[0]}=1;
}

open (IN, "KEGG_modules.tab"); ## assign kegg ID to metag module
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$anno{$a[0]}=$a[1];
	@b=split(";",$a[2]);
	for my $val (@b) {
		next unless (exists $ko{$val});
		$bme{$a[0]}{$val}=1;
	}
}
open (IN, "metab_module_assign.txt");## assign metabolite ID to metag module
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$mme{"MetaB_ME".$a[1]}{$a[0]}=1;
}
open (IN, "metab_module_description.txt"); ## metab module description
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$anno{"MetaB_ME".$a[0]}=$a[6];
}
open (IN, "ko2metabo.txt");  # internal linkage
while (<IN>) {
	chop;
	@a=split("\t",$_);
	#@b=split(";",$a[1]);
	#for my $val (@b) {
		$match{$a[0]}{$a[1]}=$a[2];
		$match2{$a[1]}{$a[0]}=$a[2];
	#}
}
open (IN, "KEGG_modules.tab"); ## KEGG metag module membership
while (<IN>) {
	chop;
	@a=split("\t",$_);
	@b=split(";",$a[2]);
	for my $val (@b) {
		$module_metag{$val}{$a[0]}=1;
	}
	$module_anno{$a[0]}=$a[1];
}
open (IN, "KEGG_cmpd2module.txt"); ## KEGG metab module membership
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$module_metab{$a[1]}{$a[0]}=1;
}
open (IN, "ko_description.txt"); ## ko description
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$desc{$a[0]}=$a[1];
}
open (IN, "metabolite_description.txt"); ##cmpd description
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$desc{$a[0]}=$a[1];
}
open (IN, "currency_cmpd.txt"); ## description for currency cmpd i.e. AMP
while (<IN>) {
	chop;
	$remove{$_}=1;
}
print "Module1\tModule2\tDescription1\tDescription2\tComponent1\tComponent2\tDesc_Comp1\tDesc_Comp2\tWays_of_Linkage\tComp_Relation\tCorrelation\tPvalue\n";
open (IN, "$sp\_metab_halla/all_associations.txt"); ## input HALLA all_associations.txt table
$dump=<IN>;
while (<IN>) {
	chop;
	@a=split("\t",$_);
	for my $key (keys %{$bme{$a[0]}}) { ## for each KO family in the metag module
		for my $key2 (keys %{$mme{$a[1]}}) { # for each metabolite in the metab module
			next if (exists $remove{$desc{$key2}});
			if (exists $match{$key}{$key2}) { ## if the correlated metabolite and KO family can be directly linked
				$total{$a[1]."\t".$a[0]."\t".$anno{$a[1]}."\t".$anno{$a[0]}."\t".$key."\t".$key2."\t".$desc{$key}."\t".$desc{$key2}."\t"."Directly_Linked"."\t".$match{$key}{$key2}."\t".$a[2]."\t".$a[3]}=1;	
			}
			else {
				for my $module (keys %{$module_metag{$key}}) {
					if (exists $module_metab{$key2}{$module}) {
						$total{$a[1]."\t".$a[0]."\t".$anno{$a[1]}."\t".$anno{$a[0]}."\t".$key."\t".$key2."\t".$desc{$key}."\t".$desc{$key2}."\t"."Direct_Same_Module"."\t".$module."\t".$a[2]."\t".$a[3]}=1;	
					}
				}
				for my $key3 (keys %{$match2{$key2}}) { # for each KO family the metabolite was linked to 
					next unless (exists $module_metag{$key3});
					for my $key4 (keys %{$module_metag{$key3}}) {
						if (exists $module_metag{$key}{$key4} and $key4 eq $a[0]) { ## if the 'linked' KO family and the 'correlated' KO family belong to the same KEGG module 
							$total{$a[1]."\t".$a[0]."\t".$anno{$a[1]}."\t".$anno{$a[0]}."\t".$key."\t".$key2."\t".$desc{$key}."\t".$desc{$key2}."\t"."Same_MetaG_Module"."\t".$key3."_".$key4."\t".$a[2]."\t".$a[3]}=1;
						}
					}	
				}
				for my $key5 (keys %{$match{$key}})  { # for each metabolite the KO family was linked to
					next unless (exists $module_metab{$key5});
					for my $key6 (keys %{$module_metab{$key5}}) {
						if (exists $module_metab{$key2}{$key6}) { ## if the 'linked' metabolite and the 'correlated' metabolite belong to the same KEGG module
							$total{$a[1]."\t".$a[0]."\t".$anno{$a[1]}."\t".$anno{$a[0]}."\t".$key."\t".$key2."\t".$desc{$key}."\t".$desc{$key2}."\t"."Same_MetaB_Module"."\t".$key5."_".$key6."\t".$a[2]."\t".$a[3]}=1;
						}
					}
				}
				for my $key7 (keys %{$match2{$key2}}) { # for each KO family the metabolite was linked to 
					for my $key8 (keys %{$match{$key}}) { # for each metabolite the KO family was linked to
						next unless (exists $module_metag{$key7} and exists $module_metab{$key8});
						for my $key9 (keys %{$module_metag{$key7}}) {
							if (exists $module_metab{$key8}{$key9}) { ## if the 'linked' KO family and the 'linked' metabolite for the 'correlated' metabolite-KO pair, belong to the same KEGG module
								$total{$a[1]."\t".$a[0]."\t".$anno{$a[1]}."\t".$anno{$a[0]}."\t".$key."\t".$key2."\t".$desc{$key}."\t".$desc{$key2}."\t"."Same_MetaG_MetaB_Module"."\t".$key7."_".$key8."_".$key9."\t".$a[2]."\t".$a[3]}=1;
							}
						}
					}
				}
			}
		}
	}
}
for my $key (keys %total) {
	print $key."\n";
}
