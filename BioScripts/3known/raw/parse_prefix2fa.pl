#####################
##get the miRBase 3-prefix of your annalysis
#special for fetch plant and annimals
#awk -F"\t" '{split($4,a,";");if(a[1]=="Metazoa"){print $1}}' /WPS/RNA/pub/database/miRBase/miRBase19/organisms.txt >/WPS/RNA/pub/database/miRBase/miRBase19/animal.prefix
#awk -F"\t" '{split($4,a,";");if(a[1]=="Viridiplantae"){print $1}}' /WPS/RNA/pub/database/miRBase/miRBase19/organisms.txt >/WPS/RNA/pub/database/miRBase/miRBase19/plant.prefix
#perl  parse_prefix2fa.pl animal.prefix mature.fa >mature_animal.fa
#perl  parse_prefix2fa.pl plant.prefix mature.fa >mature_plant.fa
##JXX 2012/08/20
use warnings;
use strict;

die "perl $0 prefix_file miRBase.fa" unless (@ARGV == 2);

my %hash=();
my $pre_file=shift;
my $infile=shift;

open(PRE,"$pre_file")||die "can't open $pre_file\n";
while(<PRE>){
	chomp;
	$hash{$_}=1;
}
close(PRE);

my $marker=0;
open(FILE,$infile);
while(<FILE>){
	chomp;
	my $tmp=$_;
	if(/^>/){
		s/\-.*//;
		s/^>//;
		if(defined($hash{$_})){
			$marker=1;
			print "$tmp\n";
		}
		else{
			$marker=0;
		}
	}
	elsif($marker){
		print "$tmp\n";
	}

}
close(FILE);
