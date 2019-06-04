#!/usr/bin/perl
use strict;
use warnings;
#-----------------------------------
# extract cdna from fa and gtf files
# 2012-09-28	
#-----------------------------------
die "perl $0 <gtf> <FA> <OUT>" unless (@ARGV==3);
open GTF, "$ARGV[0]" or die $!;
open FA, "$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]" or die $!;

my %gtf;
my %id;
while(<GTF>){
	chomp;
	my @tmp = split /\t/;
    next unless($tmp[2] =~ /exon/);
    my $gid = $1 if($tmp[8] =~ /gene_id "(.*?)";/);
    my $transID = $1 if($tmp[8] =~ /transcript_id "(.*?)";/);
    $id{$transID} = $gid;
	push @{$gtf{$tmp[0]}{$transID}},[$tmp[3], $tmp[4], $tmp[6]];
}

my %fa;
my $chr;
while(<FA>){
	chomp;
	if(/^>(\S+)/){
		$chr = $1;
		next;
	}
	$fa{$chr} .= $_;
}
sub sort_exon{
	$a->[0] <=> $b->[0]
	or $a->[1] <=> $b->[1]
}
my %geneSeq;
foreach my $key1(sort {$a cmp $b} keys %gtf){
	foreach my $key2(sort {$a cmp $b} keys %{$gtf{$key1}}){
        my @tmp=@{$gtf{$key1}{$key2}};
	my @array=sort sort_exon @tmp;
        my $seq;
        my $strand;
        for(my $i=0; $i<=$#array; $i++){
		    my $sta = $array[$i][0];
		    my $end = $array[$i][1];
		    $strand = $array[$i][2]; 
		    my $length = $end - $sta + 1;
		    if($fa{$key1}){
			if(($strand eq "+") || ($strand eq ".") || ($strand eq "-")){
				$seq .= substr ($fa{$key1}, $sta-1, $length);
        }
		   }
     }
      if($strand eq "-"){Complement_Reverse(\$seq);}
        my $geneID = $id{$key2};
        ${$geneSeq{$key1}{$geneID}} = "$seq" unless($geneSeq{$key1}{$geneID});
        my $cmpa = length($seq);
        my $cmpb = length(${$geneSeq{$key1}{$geneID}});
        ${$geneSeq{$key1}{$geneID}} = "$seq" if($cmpa >= $cmpb);
    }
}

my %tmphash;
foreach my $genekey1(sort {$a cmp $b} keys %geneSeq){
    foreach my $genekey2(sort {$a cmp $b} keys %{$geneSeq{$genekey1}}){
        my $geneSequence = "${$geneSeq{$genekey1}{$genekey2}}";
        unless(exists $tmphash{$genekey2}){
            print OUT ">$genekey2\n$geneSequence\n";
        }
        $tmphash{$genekey2}=0;
    }
}

close GTF;
close FA;
close OUT;

sub Complement_Reverse{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
		$$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
		$$seq_p=reverse($$seq_p);  
	}
}

