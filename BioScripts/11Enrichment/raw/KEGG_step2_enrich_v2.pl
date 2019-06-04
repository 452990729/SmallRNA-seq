#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use FindBin '$Bin';
use Config::Tiny;

my ($input_id, $output_dir, $species, $sample_name, $correct_method, $BlastResult);
GetOptions(
	"id=s"=>\$input_id,
	"out-dir=s"=>\$output_dir,
	"species=s"=>\$species,
	"blast-result=s" => \$BlastResult,
	"sample-names=s"=>\$sample_name,
	"correct-method=s"=>\$correct_method,
);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $R_v2153 = $Config->{srnaenv}->{R_v2153};


my $runKOBAS="$Bin/runKOBAS_v3.pl";
my $add_id="$Bin/add_id_v2.pl";
my $extract="$Bin/extract_top20.pl";
my $Rscript="$R_v2153/Rscript";
my $plot="$Bin/Pathwayscatter.R";
$correct_method||= "BH";


print "$perlExec $runKOBAS -i $input_id -o $output_dir -s $species -blast-result $BlastResult --sample-name $sample_name -f $correct_method\n";
print "sh $output_dir/annot.sh\n";
print "sh $output_dir/identify.sh\n";
print "$perlExec $add_id $output_dir/temp/$sample_name\.annotate $output_dir/$sample_name\.identify\.xls \>$output_dir/add\.$sample_name\.identify\.xls\n";
print "$perlExec $extract $output_dir/add\.$sample_name\.identify\.xls >$output_dir/top\_20\.add\.$sample_name\.identify\.xls\n";
print "sed -i \"s/'/ /g\" $output_dir/top\_20\.add\.$sample_name\.identify\.xls\n";
print "${Rscript} $plot $sample_name $output_dir $output_dir/top\_20\.add\.$sample_name\.identify\.xls\n";
