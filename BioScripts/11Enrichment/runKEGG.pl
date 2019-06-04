use strict;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
use Data::Dumper;
use FindBin '$Bin';
use Config::Tiny;


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v276};
my $kobas_v20120208 = $Config->{database}->{kobas_v20120208};
my $blast = $Config->{software}->{blast2};
my $dirDiffGeneID;
my $dirGTF;
my $eValue = 1e-10;
my $numberProcessor = 10;
my $existGenome = "n";
my $dirOutputAssemble;
my $dirScriptGetGene = "$Bin/getGene.pl";
my $dirScriptextract = "$Bin/extractcDNAfromFA.pl";
my $dirGenome;
my $dirScriptGetDiffGene = "$Bin/getDiffGene.pl";
my $dirKOBAS_generate_run = "$Bin/KOBAS_generate_run.pl";
my $dirScriptKOBAS = "$Bin/runKOBAS.pl";
my $dirSpeciesList = "$kobas_v20120208/sp.list";
my $dirSpeciesData = "$kobas_v20120208/seq_pep_v2/";
my $samplename ||="Samplename";
my $outdir=`pwd`;
$outdir =~ s/\n//;
my $dirOutputKOBAS ||=$outdir."/KOBASout";
my $blastProgram ||= "$blast/blastx";
my $abbrForSpecies ||= "ko";
my $FDRMethod ||= "BH";
my $database ||= "K/B";


GetOptions ("directory-gtf=s" => \$dirGTF,
            "gene-id=s" => \$dirDiffGeneID,
            "e=f" => \$eValue,
            "n=i" => \$numberProcessor,
            "exist-genome=s" => \$existGenome,
            "directory-assemble-result=s" => \$dirOutputAssemble,
            "directory-genome=s" => \$dirGenome,
            "directory-kobas-output=s" => \$dirOutputKOBAS,
            "directory-species-list=s" => \$dirSpeciesList,
            "directory-kegg-database=s" => \$dirSpeciesData,
			"samplename=s"=>\$samplename,
            "blast-program=s" => \$blastProgram,
            "species=s" => \$abbrForSpecies,
            "fdr=s" => \$FDRMethod,
            "db=s" => \$database,
            ) or die "Something went wrong!";

my $dirDiffGeneSeq=$dirOutputKOBAS."/"."$samplename".".diffgene.seq";
my $dirGeneSeq=$dirOutputKOBAS."/gene.seq";

unless (-d $dirOutputKOBAS){
	!system "mkdir -p $dirOutputKOBAS" or die "Something went wrong when mkdir for $dirOutputKOBAS!";
}
if ($existGenome eq "y"){
	$dirGeneSeq =~ /(.*)\//;
	unless (-d $1){
		!system "mkdir -p $1"or die "Something went wrong when mkdir for $dirGeneSeq!";
	}
#	!system "$dirScriptGetGene $dirGenome $dirGTF $dirGeneSeq" or die "Something went wrong!";
	!system "$perlExec $dirScriptextract $dirGTF $dirGenome $dirGeneSeq" or die "Something went wrong!";
	if (-s $dirGeneSeq){
		$dirDiffGeneSeq =~ /(.*)\//;
		unless (-d $1){
			!system "mkdir -p $1" or die "Something went wrong when mkdir for $dirDiffGeneSeq!";
		}
		!system "$dirScriptGetDiffGene $dirGeneSeq $dirDiffGeneID $dirDiffGeneSeq" or die "Something went wrong!";
	}
}

if ($existGenome eq "y"){
	`$pythonExec $Bin/KEGG_step1_blast.py $dirDiffGeneSeq $abbrForSpecies $dirOutputKOBAS/$samplename.xml $dirOutputKOBAS/blast.sh `;
	`sh $dirOutputKOBAS/blast.sh`;
	`$perlExec $Bin/KEGG_step2_enrich_v2.pl -id $dirDiffGeneID -out-dir $dirOutputKOBAS -species $abbrForSpecies -blast-result $dirOutputKOBAS/$samplename.xml -sample-names $samplename >$dirOutputKOBAS/run.sh`;
	`sh $dirOutputKOBAS/run.sh`;
}
elsif ($existGenome eq "n"){
	!system "$dirKOBAS_generate_run --fa $dirOutputAssemble --out-dir $dirOutputKOBAS --species $abbrForSpecies --database $database --correct-method $FDRMethod --blast-program  $blastProgram --sample-names $samplename >$dirOutputKOBAS/run.sh" or die "Something went wrong!";
	!system "sh $dirOutputKOBAS/run.sh" or die "Something went wrong!";
}
