use strict;
use warnings;
use File::Basename;
use FindBin '$Bin';
use Config::Tiny;

die "perl $0 ncRNA.unmap.fas repeat_predict_dir repeat_map_outdir" unless (@ARGV==3);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $bowtie = $Config->{software}->{bowtie1};
my $blast = $Config->{software}->{blast2};

my $small_fna=shift;	#input small-rna format file
my $in_dir=shift;        #repeat predict dir
my $out_dir=shift;	#mapping sta output dir
if(-e "$out_dir"){
	`rm -rf $out_dir`;
}

my $small_fna_fulpath=$small_fna;
if($small_fna !~ m{^/}){
	my $pwd=`pwd`;
	chomp($pwd);
	$small_fna_fulpath="$pwd/$small_fna";
}

`mkdir -p $out_dir`;
$small_fna="$out_dir/".basename($small_fna);
`ln -s $small_fna_fulpath $small_fna`;

my @types;
my %gi;

my $DR_fa=`ls $in_dir/*.DR.fna`;
chomp($DR_fa);
my $DR_gff=`ls $in_dir/*.DR.gff`;
chomp($DR_gff);
my %types_DR;
open(GFF,$DR_gff)||die "'can't open $DR_gff";
while(<GFF>){
	my @t = split("\t", $_);
	my ($id,$type) = $t[-1] =~ /^ID=([^;]+);\S+;Class=([^;]+);$/g;
	$type=~s/\//:/;
	$type=~s/\(.*\)//; #edit by jc on 2016/4/8
	$gi{$type}.=$id."\n";
#	`echo $id >>$out_dir/$type.gi`;
	$types_DR{$type}=1;
}
close(GFF);
foreach my $i(sort keys %types_DR){
	push @types,$i if $types_DR{$i}==1;
}


my $TR_fa=`ls $in_dir/*.TR.fna`;
chomp($TR_fa);
my $TR_gff=`ls $in_dir/*.TR.gff`;
chomp($TR_gff);
if(!(-e "$in_dir/DR_TR.fna")){
	`cat $DR_fa $TR_fa >$in_dir/DR_TR.fna`;
}
if(!(-e "$in_dir/DR_TR.fna.nhr")){
	`$blast/makeblastdb -in $in_dir/DR_TR.fna -dbtype nucl -parse_seqids -hash_index`;
}

my %types_TR;
open(GFF,$TR_gff)||die "'can't open $TR_gff";
while(<GFF>){
	my @t = split("\t", $_);
	my ($id,$type) = $t[-1] =~ /^ID=([^;]+);Type=([^;]+);/g;
	$type=~s/\//:/;
	$type=~s/\(.*\)//;
	$gi{$type}.=$id."\n";
#	`echo $id >>$out_dir/$type.gi`;
	$types_TR{$type}=1;
}
close(GFF);
foreach my $i(sort keys %types_TR){
	push @types,$i if $types_TR{$i}==1;
}

=head
my $IR_fa=`ls $in_dir/*.IR.fna`;
chomp($IR_fa);
if(-e $IR_fa){
	push @types,"invert_repeat";
	`cp $IR_fa $out_dir/invert_repeat.fna`;
}
=cut

my $uc_list="";
my $rc_list="";
`>$out_dir/repeat.map.fas\n`;
for my $i(0..$#types){
	if(defined($gi{$types[$i]})){
		open(OUT,">$out_dir/$types[$i].gi");
		print OUT "$gi{$types[$i]}";
		`$blast/blastdbcmd -db $in_dir/DR_TR.fna -dbtype nucl -entry_batch $out_dir/$types[$i].gi |sed \'s/lcl|//g\' >>$out_dir/$types[$i].fna`;
	}
	`echo "$bowtie/bowtie-build $out_dir/$types[$i].fna"`;
	`$bowtie/bowtie-build $out_dir/$types[$i].fna $out_dir/$types[$i].fna`;
	`$bowtie/bowtie -p 10 -v 0 -k 1 -f $out_dir/$types[$i].fna $small_fna $out_dir/$types[$i].bwt --un $out_dir/$types[$i].unmap.fas 2>$out_dir/$types[$i].mapping.log`;
	if(!(-z "$out_dir/$types[$i].bwt")){
		`$perlExec $Bin/bwt12collapse.pl $out_dir/$types[$i].bwt >>$out_dir/repeat.map.fas\n`;
		`$perlExec $Bin/genebwt12count.pl -i $out_dir/$types[$i].bwt -r $out_dir/$types[$i].fna -t $types[$i] -o $out_dir -u -s`;
		`awk \'{if(NR==1){print \"Types\\t\"\$0}else if(NR==2){print \"Mapped reference\\t\"\$0}}\' $out_dir/$types[$i].mapref.stat >$out_dir/$types[$i].map.stat`;
		`awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $out_dir/$types[$i].uc.stat >>$out_dir/$types[$i].map.stat`;
		`awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped total sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $out_dir/$types[$i].rc.stat >>$out_dir/$types[$i].map.stat`;
		$uc_list=$uc_list." $out_dir/$types[$i].uc.stat";
		$rc_list=$rc_list." $out_dir/$types[$i].rc.stat";
		$small_fna="$out_dir/$types[$i].unmap.fas";
	}

}
`cp $small_fna $out_dir/repeat.unmap.fas`;
`$perlExec $Bin/paste_col.pl $uc_list >$out_dir/uc.stat`;
`$perlExec $Bin/paste_col.pl $rc_list >$out_dir/rc.stat`;
`awk '{if(NR==1){num=NF-1;print}else if(\$1~/+\$/){for(i=1;i<=num;i++){total1[i]+=\$(i+1)}}else if(\$1~/-\$/){for(i=1;i<=num;i++){total2[i]+=\$(i+1)}}else{for(i=1;i<=num;i++){total[i]+=\$(i+1)}}}END{printf("repeat");for(i=1;i<=num;i++){printf("\\t"total[i])}printf("\\n");printf("repeat:+");for(i=1;i<=num;i++){printf("\\t"total1[i])}printf("\\n");printf("repeat:-");for(i=1;i<=num;i++){printf("\\t"total2[i])}printf("\\n");}' $out_dir/uc.stat >$out_dir/repeat.uc.stat`;
`awk '{if(NR==1){num=NF-1;print}else if(\$1~/+\$/){for(i=1;i<=num;i++){total1[i]+=\$(i+1)}}else if(\$1~/-\$/){for(i=1;i<=num;i++){total2[i]+=\$(i+1)}}else{for(i=1;i<=num;i++){total[i]+=\$(i+1)}}}END{printf("repeat");for(i=1;i<=num;i++){printf("\\t"total[i])}printf("\\n");printf("repeat:+");for(i=1;i<=num;i++){printf("\\t"total1[i])}printf("\\n");printf("repeat:-");for(i=1;i<=num;i++){printf("\\t"total2[i])}printf("\\n");}' $out_dir/rc.stat >$out_dir/repeat.rc.stat`;
my $cmd="cd $out_dir/\n$perlExec $Bin/stat_image_txt.pl -r uc.stat\n$perlExec $Bin/image_txt_bar.pl uc.stat.image.txt uc.stat Uniq\n$perlExec $Bin/stat_image_txt.pl -r rc.stat\n$perlExec $Bin/image_txt_bar.pl rc.stat.image.txt rc.stat Total\ncp repeat.* ../\ncp *.stat.FMT ../\ncp *.pdf ../\ncp *.png ../\n";
#print $cmd;
`$cmd`;
=head
`cd $out_dir/`;
`perl $Bin/stat_image_txt.pl -r uc.stat`;
`perl $Bin/image_txt_bar.pl uc.stat.image.txt repeat.uc.stat Uniq`;
`perl $Bin/stat_image_txt.pl -r rc.stat`;
`perl $Bin/image_txt_bar.pl rc.stat.image.txt repeat.rc.stat Total`;
`cp repeat.* ../`;
`cp *.stat.FMT ../`;
=cut
