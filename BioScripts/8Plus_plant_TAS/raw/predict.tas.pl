use strict;
use warnings;
use FindBin '$Bin';
use Config::Tiny;

die "perl $0 genome QC_directory" unless (@ARGV == 2);

my $genome = shift;
my $QC_dir = shift;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $srnatoolscli = $Config->{software}->{srnatoolscli};

`export PATH=$srnatoolscli/bin:\$PATH`;
if(!(-e "phase.predict")){
	`mkdir -p phase.predict`;
}
`cat $QC_dir/*/clean_data/*_clean_total.fa >phase.predict/clean_total.fa`;

`$perlExec $srnatoolscli/srna-tools.pl --tool phasing --genome $genome --srna_file phase.predict/clean_total.fa  --abundance 3 --pval 0.001 --minsize 20 --maxsize 26 --trrna --out phase.predict >phase.predict/phase.log 2>phase.predict/phase.error`;
`$perlExec $Bin/parse.tas.pl $genome phase.predict phase.out`;
