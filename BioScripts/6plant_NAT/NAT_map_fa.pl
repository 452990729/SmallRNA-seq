use strict;
use warnings;
use FindBin '$Bin';
use Config::Tiny;

sub usage{
print STDERR <<USAGE;
=========================================================================
example: perl $0 cis-NAT:cis-NAT.fa,trans-NAT:trans-NAT.fa repeat.unmap.fas sit
qsub -cwd -V -l vf=200M  sit_NAT_map.sh
=========================================================================
USAGE
exit 0;
}

&usage unless (@ARGV==3);
#die "perl $0 type1:fna1,type2:fna2,type3:fna3 reads.fa outprj\n" unless (@ARGV==3);
##type1,type2,type3,type4 as the classify order

my @input=split(",",shift);
my $read=shift;
my $prj=shift;
my $pwd = `pwd`;
$pwd =~ s/\s+//g;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $bowtie = $Config->{software}->{bowtie1};

open OUT, ">$prj\_NAT_map.sh" or die $!;

if(!(-e "$pwd/$prj")){
	print OUT "mkdir -p $pwd/$prj\n";
}

if(!(-e "$pwd/$prj/input")){
	print OUT "mkdir -p $pwd/$prj/input\n";
}

if(!(-e "$pwd/$prj/output")){
	print OUT "mkdir -p $pwd/$prj/output\n";
}

print OUT "echo Start Time:\ndate\n";
print OUT ">$pwd/$prj/output/NAT.map.fas\n";
my ($uc_list,$rc_list);
foreach my $i(0..$#input){
	my ($type,$fna)=split(":",$input[$i]);
	print OUT "echo ================== Run $type mapping ==================\n";
	print OUT "ln -sf ",get_ful_path($fna)," $pwd/$prj/input/$type.fna\n";
	print OUT "$bowtie/bowtie-build $pwd/$prj/input/$type.fna $pwd/$prj/input/$type.fna\n";
	print OUT "$bowtie/bowtie -p 8 -v 0 -k 1 -f $pwd/$prj/input/$type.fna $read $pwd/$prj/output/$type.bwt --un $pwd/$prj/output/$type.unmap.fas 2>$pwd/$prj/output/$type.mapping.log\n";
	print OUT "echo ================== Stat $type mapping  ==================\n";
	print OUT "if [[ ! -z $pwd/$prj/output/$type.bwt  ]]; then\n";
	print OUT "$perlExec $Bin/bwt12collapse.pl $pwd/$prj/output/$type.bwt >$pwd/$prj/output/$type.map.fas\n";
	print OUT "cat $pwd/$prj/output/$type.map.fas >>$pwd/$prj/output/NAT.map.fas\n";
	print OUT "$perlExec $Bin/genebwt12count.pl -i $pwd/$prj/output/$type.bwt -r $pwd/$prj/input/$type.fna -t $type -o $pwd/$prj/output/ -u -s\n";
	print OUT "awk \'{if(NR==1){print \"Types\\t\"\$0}else if(NR==2){print \"Mapped reference\\t\"\$0}}\' $pwd/$prj/output/$type.mapref.stat >$pwd/$prj/output/$type.map.stat\n";
	print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $pwd/$prj/output/$type.uc.stat >>$pwd/$prj/output/$type.map.stat\n";
	print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped total sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $pwd/$prj/output/$type.rc.stat >>$pwd/$prj/output/$type.map.stat\n";
	$uc_list.="$pwd/$prj/output/$type.uc.stat ";
	$rc_list.="$pwd/$prj/output/$type.rc.stat ";
	print OUT "fi\n";
	$read="$pwd/$prj/output/$type.unmap.fas";
}
print OUT "$perlExec $Bin/paste_col.pl $uc_list >$pwd/$prj/output/uc.stat\n";
print OUT "$perlExec $Bin/paste_col.pl $rc_list >$pwd/$prj/output/rc.stat\n";
print OUT "awk \'{if(NR==1){num=NF-1;print}else if(\$1~/+\$/){for(i=1;i<=num;i++){total1[i]+=\$(i+1)}}else if(\$1~/-\$/){for(i=1;i<=num;i++){total2[i]+=\$(i+1)}}else{for(i=1;i<=num;i++){total[i]+=\$(i+1)}}}END{printf(\"NAT\");for(i=1;i<=num;i++){printf(\"\\t\"total[i])}printf(\"\\n\");printf(\"NAT:+\");for(i=1;i<=num;i++){printf(\"\\t\"total1[i])}printf(\"\\n\");printf(\"NAT:-\");for(i=1;i<=num;i++){printf(\"\\t\"total2[i])}printf(\"\\n\");}\' $pwd/$prj/output/uc.stat >$pwd/$prj/output/NAT.uc.stat\n";
print OUT "awk \'{if(NR==1){num=NF-1;print}else if(\$1~/+\$/){for(i=1;i<=num;i++){total1[i]+=\$(i+1)}}else if(\$1~/-\$/){for(i=1;i<=num;i++){total2[i]+=\$(i+1)}}else{for(i=1;i<=num;i++){total[i]+=\$(i+1)}}}END{printf(\"NAT\");for(i=1;i<=num;i++){printf(\"\\t\"total[i])}printf(\"\\n\");printf(\"NAT:+\");for(i=1;i<=num;i++){printf(\"\\t\"total1[i])}printf(\"\\n\");printf(\"NAT:-\");for(i=1;i<=num;i++){printf(\"\\t\"total2[i])}printf(\"\\n\");}\' $pwd/$prj/output/rc.stat >$pwd/$prj/output/NAT.rc.stat\n";
print OUT "cp $read $pwd/$prj/output/NAT.unmap.fas\n";
print OUT "echo End Time:\ndate\n";

close(OUT);

sub get_ful_path{
	my $in=shift;
	if($in !~ /^\//){
		my $t=`pwd`;
		chomp($t);
		return "$t/$in";
	}else{
		return $in;
	}
}
