use strict;
use warnings;
use Getopt::Std;

# Get parameters
our %opt;
&getopts ('r:h', \%opt);
unless ($opt{r}) {
    &HELP_MESSAGE();
    exit 1;
}
if ($opt{h}) {
    &HELP_MESSAGE();
    exit 0;
}


open (RC,"$opt{r}") or die;
open (OUT1,">$opt{r}.FMT");
open (OUT2,">$opt{r}.image.txt");
##############get sample name###########
my $head=<RC>;
chomp($head);
my @heads=split("\t",$head);
my @samples;
print OUT1 "Types";
print OUT2 "Types";
foreach my $i(1..$#heads){
	print OUT1 "\t$heads[$i]";
	print OUT2 "\t$heads[$i]";
	push(@samples,$heads[$i]);
}
print OUT1 "\n";
print OUT2 "\n";

my (%hash,%hash1,%hash2);
my (%total,%total1,%total2);
my (%ambi,%ambi1,%ambi2);
my @types;
my $type="Types";
while(<RC>){
	my $ambi_marker=0;
	my @tmp=split("\t",$_);
	my $strand="";
	if($tmp[0] =~ /\?/ || $tmp[0] =~ /Unknown/){$ambi_marker=1;}
	elsif(!/^$type/){
		if(/:/){
			my @temp=split(":",$tmp[0]);
			$type=$temp[0];
		}else{$type=$tmp[0];}
		push(@types,$type);
	}
	if($tmp[0] =~ /\+$/){
		$strand="+";
	}elsif($tmp[0] =~ /\-$/){
		$strand="-";
	}
	if($strand){
		foreach my $i(0..$#samples){
			$total{$samples[$i]}+=$tmp[$i+1];
			if($strand eq "+"){$total1{$samples[$i]}+=$tmp[$i+1];}
			elsif($strand eq "-"){$total2{$samples[$i]}+=$tmp[$i+1];}
		}
		if($ambi_marker){
			foreach my $i(0..$#samples){
				$ambi{$samples[$i]}+=$tmp[$i+1];
				if($strand eq "+"){$ambi1{$samples[$i]}+=$tmp[$i+1];}
				if($strand eq "-"){$ambi2{$samples[$i]}+=$tmp[$i+1];}
			}
		}else{
			foreach my $i(0..$#samples){
				$hash{$samples[$i]}{$type}+=$tmp[$i+1];
				if($strand eq "+"){$hash1{$samples[$i]}{$type}+=$tmp[$i+1];}
				elsif($strand eq "-"){$hash2{$samples[$i]}{$type}+=$tmp[$i+1];}
			}
		}
	}
}


foreach my $i(@types){
	print OUT1 $i;
	foreach my $j(@samples){
		print OUT1 "\t",$hash{$j}{$i};
	}
	print OUT1 "\n",$i,":+";
	print OUT2 $i,":+";
	foreach my $j(@samples){
		print OUT1 "\t",$hash1{$j}{$i};
		print OUT2 "\t",$hash1{$j}{$i};
	}
	print OUT1 "\n",$i,":-";
	print OUT2 "\n",$i,":-";
	foreach my $j(@samples){
		print OUT1 "\t",$hash2{$j}{$i};
		print OUT2 "\t",$hash2{$j}{$i};
	}
	print OUT1 "\n";
	print OUT2 "\n";
}
print OUT1 "ambi";
print OUT2 "ambi";
foreach my $j(@samples){
	if(defined($ambi{$j})){
		print OUT1 "\t",$ambi{$j};
		print OUT2 "\t",$ambi{$j};
	}else{
		print OUT1 "\t0";
		print OUT2 "\t0";
	}
}
print OUT2 "\n";
close(OUT2);
print OUT1 "\nambi:+";
foreach my $j(@samples){
	if(defined($ambi1{$j})){print OUT1 "\t",$ambi1{$j};}
	else{print OUT1 "\t0";}
}
print OUT1 "\nambi:-";
foreach my $j(@samples){
	if(defined($ambi2{$j})){print OUT1 "\t",$ambi2{$j};}
	else{print OUT1 "\t0";}
}
print OUT1 "\ntotal";
foreach my $j(@samples){
        print OUT1 "\t",$total{$j};
}
print OUT1 "\ntotal:+";
foreach my $j(@samples){
        print OUT1 "\t",$total1{$j};
}
print OUT1 "\ntotal:-";
foreach my $j(@samples){
        print OUT1 "\t",$total2{$j};
}
print OUT1 "\n";
close(OUT1);

sub HELP_MESSAGE {
    my ($prog_name) = $0=~m!([^/]+)$!;
    print <<EOD;
Usage:
       perl $prog_name [args] <query file>
             -r		rc.stat/uc.stat
             -h  help
EOD
}
