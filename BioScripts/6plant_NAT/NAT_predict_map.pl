#Edit by zhangyu, 2013.09.18
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Config::Tiny;


sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description     Creat sh to run NAT-siRNA predict and mapping analysis.
Version:1.0
Usage: perl $0 [options]
       qsub -cwd -V -l vf=2g -l p=4 runNAT_predict_map.sh
Options
[mandatory parameters]
                -q  <s>     :  Unmapped fasta file from last step
                -o  <s>     :  abbr of specie
                -h|?|help   :  Show this help
[optional parameters]: choose one of the combination: -gene/-gff, -gff/-genome, -gtf/-gene, -gtf/-genome
                -gene [s]   :  Gene fasta file
		-genome[s]  :  Genome fasta file
                -gtf  [s]   :  Reference genome gtf file form Ensembl
                -gff  [s]   :  Reference genome gff file contains "gene, mRNA, CDS"
		-outdir [s] :  output directory
=========================================================================
USAGE
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
#my $blast = $Config->{software}->{blast2};

#=================================================
my ($help,$query,$abbr,$gene,$gtf,$gff,$genome,$outdir);
GetOptions(
        "h|?|help"=>\$help,
        "q=s"=>\$query,
        "o=s"=>\$abbr,
        "gene=s"=>\$gene,
        "gtf=s"=>\$gtf,
        "gff=s"=>\$gff,
	"genome=s"=>\$genome,
	"outdir=s"=>\$outdir,
);

#gff3 like this:
###gff-version 3
#scaffold_1      phytozome9_0    gene    7055662 7079200 .       +       .       ID=Si016055m.g;Name=Si016055m.g
#scaffold_1      phytozome9_0    mRNA    7055662 7079200 .       +       .       ID=PAC:19689856;Name=Si016055m;pacid=19689856;longest=1;Parent=Si016055m.g
#scaffold_1      phytozome9_0    CDS     7055662 7055825 .       +       0       ID=PAC:19689856.CDS.1;Parent=PAC:19689856;pacid=19689856
if(!defined($query) || !defined($abbr)){
        &usage;
        exit 0;
}

#chomp(my $dir = $outdir);
$outdir =~ s/\s+//g;

if(!-e "$outdir/prepare")
{
	mkdir "$outdir/prepare";
}

if(!-e "$outdir/predict")
{
        mkdir "$outdir/predict";
}

=head
open OUT, ">runNAT_predict_map.sh" or die $!;
print OUT "echo Start Time:\ndate\n";
print OUT "echo \"=================prepare data======================\"\n";
if(defined($gtf)){
	print OUT "$perlExec $Bin/gtf_to_gff.pl $gtf >$gtf.gff\n";
	$gff="$gtf.gff";
}
=cut

print "date +\"%D %T ->Start NAT predict map \" && \\\n";
print "echo \"=================prepare data======================\"\n";
if(defined($gtf)){
	my $tmp=`basename $gtf`;
	chomp $tmp;
	print "$perlExec $Bin/gtf_to_gff.pl $gtf >$outdir/prepare/$tmp.gff\n";
	$gff="$outdir/prepare/$tmp.gff";
}
if(!defined($gene)){
	my $tmp=`basename $gff`;
	chomp $tmp;
	print "$perlExec $Bin/gff2seq.pl  -g $gff -f gene -s $genome -o $outdir/prepare/$tmp.gene.fa\n";
	$gene="$outdir/prepare/$tmp.gene.fa";
}

=head
print OUT "echo \"=================cis-NAT $abbr======================\"\n";
print OUT "perl $Bin/bin/cis.pl -g $gff -od $outdir/predict -o $abbr.cis\n";
print OUT "awk \'{if(NR>1){print \$5\"\\n\"\$10}}\' $outdir/predict/$abbr.cis >$outdir/predict/$abbr.cis.id\n";
print OUT "perl $Bin/bin/fasta_id_catch_fasta.pl $outdir/predict/$abbr.cis.id $gene  $outdir/predict/$abbr.cis.fa\n";

print OUT "echo \"=================trans-NAT $abbr======================\"\n";
print OUT "perl $Bin/bin/trans.pl -q $gene -s $gene -bo $outdir/predict/blastn.out -g $gff -to $outdir/predict/$abbr.trans\n";
print OUT "awk \'{if(NR>1){print \$5\"\\n\"\$10}}\' $outdir/predict/$abbr.trans >$outdir/predict/$abbr.trans.id\n";
print OUT "perl $Bin/bin/fasta_id_catch_fasta.pl $outdir/predict/$abbr.trans.id $gene  $outdir/predict/$abbr.trans.fa\n";

print OUT "echo \"=================mapping and stat======================\"\n";
print OUT "perl $Bin/bin/NAT_map_fa.pl cis-NAT:$outdir/predict/$abbr.cis.fa,trans-NAT:$outdir/predict/$abbr.trans.fa $query $abbr\n";
print OUT "sh $abbr\_NAT_map.sh\n";

print OUT "echo End Time:\ndate\n";
=cut

print  "echo \"=================cis-NAT $abbr======================\"\n";
print  "$perlExec $Bin/cis.pl -g $gff -od $outdir/predict -o $abbr.cis\n";
print  "awk \'{if(NR>1){print \$5\"\\n\"\$10}}\' $outdir/predict/$abbr.cis >$outdir/predict/$abbr.cis.id\n";
print  "$perlExec $Bin/fasta_id_catch_fasta.pl $outdir/predict/$abbr.cis.id $gene  $outdir/predict/$abbr.cis.fa\n";

print  "echo \"=================trans-NAT $abbr======================\"\n";
print  "$perlExec $Bin/trans.pl -q $gene -s $gene -bo $outdir/predict/blastn.out -g $gff -to $outdir/predict/$abbr.trans\n";
print  "awk \'{if(NR>1){print \$5\"\\n\"\$10}}\' $outdir/predict/$abbr.trans >$outdir/predict/$abbr.trans.id\n";
print  "$perlExec $Bin/fasta_id_catch_fasta.pl $outdir/predict/$abbr.trans.id $gene  $outdir/predict/$abbr.trans.fa\n";

print  "echo \"=================mapping and stat======================\"\n";
print  "$perlExec $Bin/NAT_map_fa.pl cis-NAT:$outdir/predict/$abbr.cis.fa,trans-NAT:$outdir/predict/$abbr.trans.fa $query $abbr\n";
print  "sh $abbr\_NAT_map.sh\n";

#print OUT "echo End Time:\ndate\n";
print "date +\"%D %T ->Finish NAT predict map \"\n";
sub get_ful_path{
        my $in=shift;
        if($in !~ /^\//){
                my $t=`pwd`;
                chomp($t);
                return "$t/$in";
        }
        else{
                return $in;
        }
}
