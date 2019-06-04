use strict;
use warnings;
use FindBin '$Bin';
use Config::Tiny;

sub usage{
print STDERR <<USAGE;
=========================================================================
example: perl $0 exon:exon.fa,intron:intron.fa former.unmap.fas ssc
qsub ssc_gene_map.sh

format:format.fa can more composition
=========================================================================
USAGE
exit 0;
}

die "perl $0 type1:fna1,type2:fna2,type3:fna3 reads.fa outprj\n" unless (@ARGV==3);
##type1,type2,type3,type4 as the classify order

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $bowtie1 = $Config->{software}->{bowtie1};
my $bowtie2 = $Config->{software}->{bowtie2};

my @input=split(",",shift);
my $read=shift;
my $prj=shift;
my $bowtie_soft="";
my $bowtie_build_soft="";
my $pwd = `pwd`;
$pwd =~ s/\s+//g;

open OUT, ">$prj\_gene_map.sh" or die $!;

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
print OUT ">$pwd/$prj/output/gene.map.fas\n";
my ($uc_list,$rc_list);
foreach my $i(0..$#input){
	my ($type,$fna)=split(":",$input[$i]);
        my $fna_path=get_ful_path($fna);
	my $a=`du -s $fna_path`;
	my $size=(split("\t",$a))[0];
	if( $size >2000000 ){
        	$bowtie_soft = "$bowtie2/bowtie --large-index";
        	$bowtie_build_soft = "$bowtie2/bowtie-build --large-index";
        }else{
		$bowtie_soft="$bowtie1/bowtie";
		$bowtie_build_soft="$bowtie1/bowtie-build";
} 
	print OUT "echo ================== Run $type mapping ==================\n";
	print OUT "ln -sf ",get_ful_path($fna)," $pwd/$prj/input/$type.fna\n";
	print OUT "$bowtie_build_soft $pwd/$prj/input/$type.fna $pwd/$prj/input/$type.fna\n";
	print OUT "$bowtie_soft -p 10 -v 0 -k 1 -f $pwd/$prj/input/$type.fna $read $pwd/$prj/output/$type.bwt --un $pwd/$prj/output/$type.unmap.fas 2>$pwd/$prj/output/$type.mapping.log\n";
	print OUT "echo ================== Stat $type mapping  ==================\n";
	print OUT "if [[ ! -z $pwd/$prj/output/$type.bwt  ]]; then\n";
	print OUT "$perlExec $Bin/bwt12collapse.pl $pwd/$prj/output/$type.bwt >$pwd/$prj/output/$type.map.fas\n";
	print OUT "cat $pwd/$prj/output/$type.map.fas >>$pwd/$prj/output/gene.map.fas\n";
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
print OUT "echo End Time:\ndate\n";
print OUT "cp $read $pwd/$prj/output/gene.unmap.fas\n";
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
