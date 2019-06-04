#use strict;
#use warnings;

die "perl $0 uc/rc.stat......" unless (@ARGV >0);

my @samples;
my @rownames;
my %hash;
foreach my $file(@ARGV){
	open(IN,"$file") ;
	my $header=<IN>;
	chomp($header);
	my @heads=split("\t",$header);
	if(@samples==0){
		foreach my $i(1..$#heads){
        		push(@samples,$heads[$i]);
		}
	}else{
		foreach my $i(1..$#heads){
			my $marker=1;
			foreach my $j(0..$#samples){
				if($heads[$i] eq $samples[$j]){$marker=0;last;}
			}
			if($marker){push (@samples,$heads[$i]);}
		}
	}	
#print join("\t",@samples),"\n";
	while(<IN>){
		chomp;
		my @tmp=split("\t",$_);
		push(@rownames,$tmp[0]);
		foreach my $i(1..$#tmp){
			$hash{$tmp[0]}{$heads[$i]}=$tmp[$i];
		}
	}
	close(IN);

}

print "Types";
foreach my $i(sort { $a cmp $b} @samples){
	print "\t$i";
#	print "\t$samples[$i]";
}
print "\n";

foreach my $i(@rownames){
	print "$i";
#	print "$rownames[$i]";
	foreach my $j(sort { $a cmp $b} @samples){
		if(defined($hash{$i}{$j})){print "\t$hash{$i}{$j}";}
		else{print "\t0";}
#		print "\t$hash{$rownames[$i]}{$samples[$j]}";
	}
	print "\n";
}
