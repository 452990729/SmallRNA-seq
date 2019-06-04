#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

my $input_fa;
my $output_dir;
my $species;
my $sample_name;
my $database;
my $correct_method;
my $blast_program;
GetOptions(
"fa=s"=>\$input_fa,
"out-dir=s"=>\$output_dir,
"species=s"=>\$species,
"sample-names=s"=>\$sample_name,
"database=s"=>\$database,
"correct-method=s"=>\$correct_method,
"blast-program=s"=>\$blast_program);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $R_v2153 = $Config->{srnaenv}->{R_v2153};
my $runKOBAS="$Bin/runKOBAS.pl";
my $add_id="$Bin/add_id.pl";
my $extract="$Bin/extract_top20.pl";
my $Rscript="$R_2153/Rscript";
my $plot="$Bin/Pathwayscatter.R";

print "$perlExec $runKOBAS -i $input_fa -o $output_dir -s $species --sample-name $sample_name -d $database -f $correct_method -bp $blast_program\n";
print "sh $output_dir/blast.sh\n";
print "sh $output_dir/annot.sh\n";
print "sh $output_dir/identify.sh\n";
print "$perlExec $add_id $output_dir/temp/$sample_name\.annotate $output_dir/$sample_name\.indentify\.xls \>$output_dir/add\.$sample_name\.indentify\.xls\n";
print "$perlExec $extract $output_dir/add\.$sample_name\.indentify\.xls >$output_dir/top\_20\.add\.$sample_name\.indentify\.xls\n";
print "sed -i \"s/'/ /g\" $output_dir/top\_20\.add\.$sample_name\.indentify\.xls\n";
print "${Rscript} $plot $sample_name $output_dir $output_dir/top\_20\.add\.$sample_name\.indentify\.xls\n";
