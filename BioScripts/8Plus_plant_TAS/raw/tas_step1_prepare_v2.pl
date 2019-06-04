use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long;
use File::Basename;
use Config::Tiny;

sub usage{
print STDERR <<USAGE;
=========================================================================
example: perl $0 -genome genome.fa -dir QC_dir -out prj
qsub -cwd -V -l vf=200M prj_tas_step1.sh
=========================================================================
options:
	-h|?|--help		help information
	-genome			genome sequence(fasta)
	-dir			QC directory
	-out			project name
=========================================================================
USAGE
exit 0;
}

my($help,$genome,$QC_dir,$prj);
GetOptions(
	"h|?|help"=>\$help,
	"genome=s"=>\$genome,
	"dir=s"=>\$QC_dir,
	"out=s"=>\$prj,
);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $TAS_db = $Config->{database}->{known_TAS};
my $srnatoolscli = $Config->{software}->{srnatoolscli};

if(!defined($genome) || !defined($QC_dir) || !defined($prj) || defined($help)){
	&usage;
}

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$genome=get_ful_path($genome);
$QC_dir=get_ful_path($QC_dir);
$prj=get_ful_path($prj);

open OUT, ">$prj\_tas_step1.sh" or die $!;
print OUT "echo Start Time:\ndate\n";
if(!(-e "$prj")){
	print OUT "mkdir -p $prj\n";
}
print OUT "echo ==================Fetch known_TAS ==================\n";
print OUT "cd $prj\n";
print OUT "$perlExec $Bin/known.tas_v2.pl $genome $TAS_db\n";
print OUT "echo ==================Fetch novel_TAS ==================\n";
print OUT "cd $prj\n";
print OUT "export PATH=$srnatoolscli/bin:\$PATH\n";
if(-e "$genome.nhd"){
	print OUT "$perlExec $Bin/predict.tas.pl $genome $QC_dir\n";
}
else{
	my $file=basename($genome);
	print OUT "$perlExec $Bin/predict.tas.pl $prj/known_TAS/ref/$file $QC_dir\n";
}
print OUT "echo ==================CAT known_TAS novel_TAS ==================\n";
print OUT "cd $prj\n";
print OUT "cat $prj/known_TAS/known_TAS.blastn.fa $prj/phase.out/novel_TAS.phased.fa >$prj/TAS.phased.fa\n";
print OUT "cd $pwd\n";
print OUT "echo End Time:\ndate\n";

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

