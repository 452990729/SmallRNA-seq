use strict;
use Getopt::Std;

## options
my %options=();
getopts("r:i:t:o:usW",\%options);

my $usage="usage:
\tperl $0 [options] -i read.ref.bwt -r ref.fa -t types -o outdir
[options]
[mandatory parameters]
\t-i	read.ref.bwt
\t-r	ref.fa
\t-t	types
\t-o	outdir
[optional parameters]
\t-W    read counts are weighted by their number of mappings. e.g. A read maps twice so each position gets 0.5 added to its read profile
\t-u	output uniq statistic
\t-s	output sense/antisense statistic
\n";

if(not $options{'i'} or not $options{'r'} or not $options{'t'} or not $options{'o'}){
	die $usage;
}

my %hash;
my $type=$options{'t'};
my $out=$options{'o'};

open IN,"<$options{'r'}" or die "File $options{'r'} not found\n"; # ref.fa
while(<IN>){
	chomp;
	if(/^>(\S+)/){
		my $id = $1;  
		$hash{$id}{'read'}=0;
		$hash{$id}{'uniq'}=0;
	}
}
close(IN);

my %mapcounts;
if($options{'W'}){
	open IN,"<$options{'i'}" or die "File $options{'i'} not found\n";  # read.ref.bwt
	while(<IN>){
	#R4N_290_x1	+	ath-MIR403	102	TTAGATTCACGCACAAACTCG	IIIIIIIIIIIIIIIIIIIII	0	
		if(/^(\S+)/){$mapcounts{$1}++;}
	}
	close(IN);
}


open IN,"<$options{'i'}" or die "File $options{'i'} not found\n";
my %total_sample;
my %hash_sample;
my $sample;
my @scores;
my $len_sc;
my $uniq=1;

#if select option -s #output sense/antisense statistic  
#this 6 variations will be used
my %hash_se;		#recording each gene total sense count
my %hash_an;		#recording each gene total annti-sense count
my %hash_sample_se;	#recording each sample each gene sense count
my %hash_sample_an;	#recording each sample each gene annti-sense count
my %total_sample_se;	#recording each sample total sense count
my %total_sample_an;	#recording each sample total annti-sense count
my %ge_sample;		#recording each sample mapped gene number

while(<IN>){
	chomp;
	#R4N_290_x1	+	ath-MIR403	102	TTAGATTCACGCACAAACTCG	IIIIIIIIIIIIIIIIIIIII	0	
	my @line = split(/\t/);
	@scores = split(/_x/,$line[0]);
	$sample = $1 if($scores[0] =~ /^(\S+)_/);
	$len_sc = $scores[$#scores];
	if($options{'W'}){$len_sc=$len_sc/$mapcounts{$line[0]};$uniq=1/$mapcounts{$line[0]};}
	$hash{$line[2]}{'uniq'}+= $uniq;
	$hash{$line[2]}{'read'}+= $len_sc;
	$hash_sample{$sample}{$line[2]}{'uniq'}+= $uniq;
	$hash_sample{$sample}{$line[2]}{'read'}+= $len_sc;
	$total_sample{$sample}{'uniq'}+= $uniq;
	$total_sample{$sample}{'read'}+= $len_sc;
	if($options{'s'}){
		if($line[1] eq "+"){
			$hash_se{$line[2]}{'uniq'}+= $uniq;
			$hash_se{$line[2]}{'read'}+= $len_sc;
			$hash_sample_se{$sample}{$line[2]}{'uniq'}+= $uniq;
			$hash_sample_se{$sample}{$line[2]}{'read'}+= $len_sc;
			$total_sample_se{$sample}{'uniq'}+= $uniq;
			$total_sample_se{$sample}{'read'}+= $len_sc;
		}
		else{
			$hash_an{$line[2]}{'uniq'}+= $uniq;
			$hash_an{$line[2]}{'read'}+= $len_sc;
			$hash_sample_an{$sample}{$line[2]}{'uniq'}+= $uniq;
			$hash_sample_an{$sample}{$line[2]}{'read'}+= $len_sc;
			$total_sample_an{$sample}{'uniq'}+= $uniq;
			$total_sample_an{$sample}{'read'}+= $len_sc;
		}
	}
}
close(IN);

open OUTG,">$out/$type.rc";
print OUTG "#title\ttotal";
my $total_t=1000000;

for my $sample(sort keys %hash_sample){
	print OUTG "\t$sample";
	if($options{'s'}){print OUTG "\t$sample(+)\t$sample(-)";}
}

for my $sample(sort keys %hash_sample){
	print OUTG "\t$sample(norm)";
	if($options{'s'}){print OUTG "\t$sample(+|norm)\t$sample(-|norm)";}
}

print OUTG "\n";

for(sort keys %hash){
	print OUTG $_,"\t",$hash{$_}{'read'};
	my $marker=0;
	for my $sample(sort keys %hash_sample){
		if($hash_sample{$sample}{$_}{'read'} > 0){
			$ge_sample{$sample}++;
			$marker=1;
			print OUTG "\t",sprintf("%.2f",$hash_sample{$sample}{$_}{'read'});
		}else{
			print OUTG "\t0";
		}
		if($options{'s'}){
			if($hash_sample_se{$sample}{$_}{'read'} > 0){
				print OUTG "\t",sprintf("%.2f",$hash_sample_se{$sample}{$_}{'read'});
			}else{print OUTG "\t0";}
			if($hash_sample_an{$sample}{$_}{'read'} > 0){
				print OUTG "\t",sprintf("%.2f",$hash_sample_an{$sample}{$_}{'read'});
			}else{print OUTG "\t0";}
		}
	}
	if($marker){$ge_sample{'Total'}++;}
	for my $sample(sort keys %hash_sample){
		if($hash_sample{$sample}{$_}{'read'} > 0){
			print OUTG "\t",sprintf("%.2f",$total_t*$hash_sample{$sample}{$_}{'read'}/$total_sample{$sample}{'read'});
		}else{print OUTG "\t0";}
		if($options{'s'}){
			if($hash_sample_se{$sample}{$_}{'read'} > 0){
				print OUTG "\t",sprintf("%.2f",$total_t*$hash_sample_se{$sample}{$_}{'read'}/$total_sample_se{$sample}{'read'});
			}else{print OUTG "\t0";}
			if($hash_sample_an{$sample}{$_}{'read'} > 0){
				print OUTG "\t",sprintf("%.2f",$total_t*$hash_sample_an{$sample}{$_}{'read'}/$total_sample_an{$sample}{'read'});
			}else{print OUTG "\t0";}
		}
	}
	print OUTG "\n";
}

print OUTG "\n";
close(OUTG);

open TOTAL,">$out/$type.mapref.stat";
my @join_sample_ge;
print TOTAL "Total";
if($ge_sample{'Total'}>0){push @join_sample_ge,$ge_sample{'Total'};}
else{push @join_sample_ge,0;}
for my $sample(sort keys %hash_sample){
	print TOTAL "\t$sample";
	if($ge_sample{$sample}>0){push @join_sample_ge,$ge_sample{$sample};}
	else{push @join_sample_ge,0;}
}
print TOTAL "\n",join("\t",@join_sample_ge),"\n";
close(TOTAL);

open TOTAL,">$out/$type.rc.stat";
my @join_sample_read;
my @join_total_read;
for my $sample(sort keys %hash_sample){
	push @join_sample_read,$sample;
	if($total_sample{$sample}{'read'}>0){push @join_total_read,$total_sample{$sample}{'read'};}
	else{push @join_total_read,0;}
}
print TOTAL "Types\t",join("\t",@join_sample_read),"\n";
print TOTAL "$type\t",join("\t",@join_total_read),"\n";

if($options{'s'}){
	my @join_total_read_se;
	for my $sample(sort keys %hash_sample){
		if($total_sample_se{$sample}{'read'}>0){push @join_total_read_se,$total_sample_se{$sample}{'read'};}
		else{push @join_total_read_se,0;}
	}
	print TOTAL "$type:+\t",join("\t",@join_total_read_se),"\n";

	my @join_total_read_an;
	for my $sample(sort keys %hash_sample){
		if($total_sample_an{$sample}{'read'}>0){push @join_total_read_an,$total_sample_an{$sample}{'read'};}
		else{push @join_total_read_an,0;}
	}
	print TOTAL "$type:-\t",join("\t",@join_total_read_an),"\n";
}
close(TOTAL);

if($options{'u'}){
	open OUTG,">$out/$type.uc";
	print OUTG "#title\ttotal";
	for my $sample(sort keys %hash_sample){
		print OUTG "\t$sample";
		if($options{'s'}){print OUTG "\t$sample(+)\t$sample(-)";}
	}
	print OUTG "\n";

	for(sort keys %hash){
		print OUTG $_,"\t",$hash{$_}{'uniq'};
		for my $sample(sort keys %hash_sample){
			if($hash_sample{$sample}{$_}{'uniq'} > 0){
				print OUTG "\t",sprintf("%.2f",$hash_sample{$sample}{$_}{'uniq'});
			}else{print OUTG "\t0";}
			if($options{'s'}){
				if($hash_sample_se{$sample}{$_}{'uniq'} > 0){
					print OUTG "\t",sprintf("%.2f",$hash_sample_se{$sample}{$_}{'uniq'});
				}else{print OUTG "\t0";}
				if($hash_sample_an{$sample}{$_}{'uniq'} > 0){
					print OUTG "\t",sprintf("%.2f",$hash_sample_an{$sample}{$_}{'uniq'});
				}else{print OUTG "\t0";}
			}
		}
		print OUTG "\n";
	}

	print OUTG "\n";
	close(OUTG);

	open TOTAL,">$out/$type.uc.stat";
	my @join_sample_uniq;
	my @join_total_uniq;
	for my $sample(sort keys %hash_sample){
		push @join_sample_uniq,$sample;
		if($total_sample{$sample}{'uniq'}>0){push @join_total_uniq,$total_sample{$sample}{'uniq'};}
		else{push @join_total_uniq,0}
	}
	print TOTAL "Types\t",join("\t",@join_sample_uniq),"\n";
	print TOTAL "$type\t",join("\t",@join_total_uniq),"\n";
	
	if($options{'s'}){
		my @join_total_uniq_se;
		for my $sample(sort keys %hash_sample){
			if($total_sample_se{$sample}{'uniq'}>0){push @join_total_uniq_se,$total_sample_se{$sample}{'uniq'};}
			else{push @join_total_uniq_se,0;}
		}
		print TOTAL "$type:+\t",join("\t",@join_total_uniq_se),"\n";


		my @join_total_uniq_an;
		for my $sample(sort keys %hash_sample){
			if($total_sample_an{$sample}{'uniq'}>0){push @join_total_uniq_an,$total_sample_an{$sample}{'uniq'};}
			else{push @join_total_uniq_an,0}
		}
		print TOTAL "$type:-\t",join("\t",@join_total_uniq_an),"\n";
	}
	close(TOTAL);
}
