use strict;
use warnings;

sub usage{
print STDERR <<USAGE;
=========================================================================
example: perl $0 known/mature.readcount,novel/mature.readcount outdir/out_prefix group groupname

readcount file format:
*****************************
type\tsample1\tsample2\tsample3...\n
id\t1\t2\t3...\n
*****************************
=========================================================================
USAGE
exit 0;
}

&usage unless (@ARGV==4);
my @infile=split(",",shift);
my $out_prefix=shift;
my $group=shift;
my $groupname=shift;

my (@sample,$type);
my %count_type_sample;
my %total_count;
foreach my $i (@infile){
	open(IN,"$i");
	my $header=<IN>;
	chomp($header);
	my @tmp_header=split("\t",$header);
	if( !($type) ) {
		($type,@sample)=split("\t",$header);
	}
	while(<IN>){
		chomp;
		my @tmp=split /\t/ ;
		foreach my $j(1..$#tmp){
			$count_type_sample{$tmp[0]}{$tmp_header[$j]}=$tmp[$j];
			$total_count{$tmp_header[$j]}+=$tmp[$j];
		}
	}
	close(IN);
}

open(OUT1,">$out_prefix.readcount");
open(OUT2,">$out_prefix.tpm");
print OUT1 "sRNA\t",join("\t",@sample),"\n";
print OUT2 "sRNA\t",join("\t",@sample),"\n";
#print $type,"\t",join("\t",@sample),"\n";
foreach my $i(sort {$a cmp $b} keys %count_type_sample){
	print OUT1 $i;
	print OUT2 $i;
	foreach my $j(@sample){
		print OUT1 "\t",$count_type_sample{$i}{$j};
		print OUT2 "\t",(1000000*$count_type_sample{$i}{$j})/$total_count{$j};
	}
	print OUT1 "\n";
	print OUT2 "\n";
}



my @group_compose=split /,/,$group;
my @group_name=split /,/,$groupname;
my %group_compose_detail;
foreach my $i(0..$#group_compose){
    my @tmp=split(":",$group_compose[$i]);
    foreach my $j(@tmp){
        $group_compose_detail{$group_name[$i]}{$j}=1;
    }
}

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
my $dir = "$pwd/diffAnalysisResult";

open(IN1,"$dir/merged.readcount");
my $header1=<IN1>;
chomp($header1);
my @header1=(split/\t/,$header1); shift @header1; # 2013-10-22 for check yangzie
my %count;
my %count_sample;  # 2013-10-22 for check yangzie
my ($types1,@sample_title1)=split("\t",$header1);
while(<IN1>){
	chomp;
	my @tmp=split /\t/;
	for my $i(1..$#tmp){
		$count{$tmp[0]}{$sample_title1[$i-1]}=$tmp[$i];
	}
	$count_sample{$tmp[0]}=join("\t",@tmp[1..$#tmp]); # 2013-10-22 for check yangzie
}
close(IN1);
open(IN2,"$dir/merged.tpm");
my $header2=<IN2>;
chomp($header2);
my @header2=(split/\t/,$header2); shift @header2; # 2013-10-22 for check yangzie
my %tpm;  
my %tpm_sample; #  2013-10-22 for check yangzie
my ($types2,@sample_title2)=split("\t",$header2);
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	for my $i(1..$#tmp){
		$tpm{$tmp[0]}{$sample_title2[$i-1]}=$tmp[$i];
	}
	$tpm_sample{$tmp[0]}=join("\t",@tmp[1..$#tmp]); # 2013-10-22 for check yangzie
}
close(IN2);

foreach my $i(keys %group_compose_detail){
	my $number=scalar(keys $group_compose_detail{$i});
	if($number!=1){
		foreach my $j(keys $group_compose_detail{$i}){
			foreach my $m(keys %count){
				$count{$m}{$i}+=$count{$m}{$j}/$number;
				$tpm{$m}{$i}+=$tpm{$m}{$j}/$number;
			}
		}
	}
}
open(OUT1,">$dir/meanscount.txt");
open(OUT2,">$dir/meanstpm.txt");
=head
print OUT1 "sRNA\t",join("\t",@group_name),"\t",join("\t",@header1),"\n"; # 2013-10-22 for check yangzie
print OUT2 "sRNA\t",join("\t",@group_name),"\t",join("\t",@header2),"\n"; # 2013-10-22 for check yangzie
foreach my $i(sort {$a cmp $b} keys %count){
	print OUT1 "$i";
	print OUT2 "$i";
	foreach my $j(@group_name){
		print OUT1 "\t",$count{$i}{$j};
		print OUT2 "\t",$tpm{$i}{$j};
	}
	print OUT1 "\t$count_sample{$i}\n"; # 2013-10-22 for check yangzie
	print OUT2 "\t$tpm_sample{$i}\n"; # 2013-10-22 for check yangzie
}
close(OUT1);
close(OUT2);
=cut 
my @name;   # No need %count_sample!  2013-10-24 modify
for (sort keys %count)
{
    @name=sort keys $count{$_};
}
print OUT1 "sRNA\t",join("\t",@name),"\n";
print OUT2 "sRNA\t",join("\t",@name),"\n";
for my $sRNA(sort keys %count)
{
    print OUT1 "$sRNA";
    print OUT2 "$sRNA";
    for my $mem(@name)
    {   
        print OUT1 "\t$count{$sRNA}{$mem}";
        print OUT2 "\t$tpm{$sRNA}{$mem}";
    }   
    print OUT1 "\n";
    print OUT2 "\n";
}

