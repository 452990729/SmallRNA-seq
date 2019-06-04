#!/usr/bin/perl
die "perl $0 <ko annotation from annotate.py> <enrichment result from identify.py> > <output>\n" unless @ARGV==2;
my $annot_result=shift;
my $path_result=shift;
my %ko;
my %Entrez;

open AN, $annot_result or die "cannot open $annot_result !";
#print "open the annot file the first time\n";
while (<AN>) {
	chomp;
	if (/\/\/\/\//) {
#		print "finished the first time reading annot file\n";
		last;
	}
	if (/^\w/) {
		my @tmp=split /\t/;
		my @IDs=split /\|/,$tmp[1];
		if ($IDs[0] eq "None") {
			$ko{$tmp[0]}=" ";
			$Entrez{$tmp[0]}=" ";
		}
		else
		{
			$ko{$tmp[0]}=$IDs[0];
			$Entrez{$tmp[0]}=" ";
		}
	}
}

$/="\/\/\/\/\n";
#print "change the line break symbol\n";
while (<AN>) {
	chomp;
	my $line=$_;
	my @tmp=split /\n/,$line;
	my @tmp1=split /\s+/,$tmp[0];
	if($tmp1[0]=~/Query:/){
		if($line=~/Entrez Gene ID:\s+(?<entrez>\d+)/){
			$Entrez{$tmp1[1]}=$+{entrez};
		}
	}
	else
		{next;}
}
close AN;

$/="\n";
#print "change the line break symbol back!\n";
open PA, $path_result or die "cannot open $path_result !";
while (<PA>) {
	chomp;
	if (/^\#Term/) {
#		print "Oh, find the title\n";
		my @tmp=split /\t/;
		print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\tKEGG_ID/KO\tEntrez_ID\t$tmp[8]\n";
	}
	elsif (/^[\w@]/) {
		my @tmp=split /\t/;
		my @id=split/\|/,$tmp[7];
		my $ko_id="";
		my $Entrez_id="";
		foreach my $i (@id) {
			$ko_id .=$ko{$i}."\|";
			$Entrez_id .=$Entrez{$i}."\|";
		}
		print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$ko_id\t$Entrez_id\t$tmp[8]\n";
	}
	else
	{
		print "$_\n";
	}
}
close PA;
