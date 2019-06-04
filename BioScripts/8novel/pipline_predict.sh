#!/bin/bash

shopt -s -o nounset


function USAGE {
	echo ""
	echo "Usage: $0 -o prefix -r reference.fasta [options]"
	echo "    Options for miRNA prediction:"
	echo "        -o  <str>     abbreviation for project name, 3 letter code for the sequencing library or the species of interest, required"
	echo "        -i  <str>   sequence reads file in FASTA format, uniqe merged, required"
	echo "        -r  <str>     The prefix of the bowtie index, for the genome reference file"
	echo "                      For instance, if the reference is 'database.fasta', then the prefix is 'database' and building-command is:"
	echo "                      'bowtie-build -f database.fasta database'"
	echo "        -M  <str>     miRNAs mature reference, miRBase miRNA sequences in fasta format. These should be the known mature sequences"
	echo "                      for the species being analyzed."
	echo "        -s  <str>     miRBase miRNA sequences in fasta format. These should be the pooled known mature sequences for 1-5 species"
	echo "                      closely related to the species being analyzed."
	echo "        -b  <int>     minimum score cut-off for predicted novel miRNAs to be displayed, default=1"
	echo "        -c            disable randfold analysis"
	echo "        -g  <int>     maximum number of precursors to analyze when automatic excision gearing is used."
	echo "                      default=50000, if set to -1 all precursors will be analyzed"
	echo "        -u  <str>     species being analyzed - this is used to link to the appropriate UCSC browser entry"
	echo "        -t  <int>     temperature cut-off for RNAfold when calculating secondary structures of RNAs, default=22"
	echo "        -m  <1/2/3>   predicttion model, 1: animal; 2: monocot; 3 dicots.  default=1"
	echo "        -w  <str>     abbreviation of species from miRBase used to adjust bowtie version"
	echo ""
	echo "    Options for Bowtie:"
	echo "        -p  <int>   number of processors to use, default=1"
	echo "        -k  <int>   a read is allowed to map up to this number of positions in the genome, default=5; for plant, 15 is recommaned"
}

if [ $# -eq 0 ]; then
	USAGE;
	exit 192;
fi

declare -rx SCRIPT=${0##*/}
declare -r OPTSTRING="o:i:r:t:b:s:M:cg:u:m:w:v:p:n:k:h"
declare SWITCH
declare genome 
declare species=""
declare known_mature
declare other_mature=" "
declare prj
declare PROJNAME 
declare INPUT
declare abbr
declare -i TEMP=22
declare -i DROSHA=0
declare -i MIS=0
declare -i CPU=1
declare -i MAXREP=5
declare -i mode=1
declare -i score=1
declare -i rand=1
declare -i gear=50000

readIni() {
    file=$1;section=$2;item=$3;
    val=`awk -F ' = ' '/\['${section}'\]/{a=1} (a==1 && "'${item}'"==$1){a=0;print $2}' $file`
    echo ${val}
}

cpath=`dirname $0`
perlExec=`readIni $cpath/../../BioModule/Settings/config.ini srnaenv perl_v5182`
bowtie=`readIni $cpath/../../BioModule/Settings/config.ini software bowtie1`
bowtie122=`readIni $cpath/../../BioModule/Settings/config.ini software bowtie122`
MIREVO=`readIni $cpath/../../BioModule/Settings/config.ini software miREvo`
RNAfold=`readIni $cpath/../../BioModule/Settings/config.ini software RNAfold`
randfold=`readIni $cpath/../../BioModule/Settings/config.ini software randfold`

program_dir="$MIREVO/script"

while getopts "$OPTSTRING" SWITCH ; do
	case $SWITCH in
		h)  USAGE;
			exit 192;
		;;
		r) genome="$OPTARG"
		;;
		i) INPUT="$OPTARG"
		;;
		o) prj="$OPTARG"
		;;
		t) TEMP="$OPTARG"
		;;
		b) score="$OPTARG"
		;;
		M) known_mature="$OPTARG"
		;;
		s) other_mature="$OPTARG"
		;;
		c) rand=0
		;;
		g) gear="$OPTARG"
		;;
		u) species="$OPTARG"
		;;
		m) mode="$OPTARG"
		;;
		w) abbr="$OPTARG"
		;;
		v) MIS="$OPTARG"
		;;
		p) CPU="$OPTARG"
		;;
		k) MAXREP="$OPTARG"
		;;
		\?) exit 192
		;;
		*) printf "miREvo predict: $LINENO: %s\n" "script error: unhandled argument"
		exit 192
		;;
	esac
done

PROJNAME=`basename $prj`
NOMAP=$prj/predict.nomap.fas
MAP=$prj/predict.map.fas
cmdlog=$prj/predict.cmd

if [[ ! -e $prj ]]; then
	mkdir -p $prj
fi

if [[ $known_mature == " " ]]; then
	echo "Please povide a fasta file for mature sequence of your species with options -M"
	echo "Exit now"
	exit 192
fi
known_mature_new=$PROJNAME.`basename $known_mature`
echo "$perlExec $program_dir/rna2dna.pl $known_mature > $prj/$known_mature_new" > $cmdlog
$perlExec $program_dir/rna2dna.pl $known_mature > $prj/$known_mature_new

if [[ $other_mature != " " ]]; then
other_mature_new=$PROJNAME.`basename $other_mature`
echo "$perlExec $program_dir/rna2dna.pl $other_mature > $prj/$other_mature_new" >> $cmdlog
$perlExec $program_dir/rna2dna.pl $other_mature > $prj/$other_mature_new
fi

echo "Bowtie mapping ..." # 4 lines
if [[ $abbr =~ "tae" ]];then
   echo "$bowtie122/bowtie --large-index -p $CPU -f -n 0 -e 80 -l 18 -a -m $MAXREP --best --strata $genome $INPUT > $prj/predict.mapping.bwt" >> $cmdlog
   $bowtie122/bowtie --large-index -p $CPU -f -n 0 -e 80 -l 18 -a -m $MAXREP --best --strata $genome $INPUT > $prj/predict.mapping.bwt
else
   echo "$bowtie/bowtie -p $CPU -f -n 0 -e 80 -l 18 -a -m $MAXREP --best --strata $genome $INPUT > $prj/predict.mapping.bwt" >> $cmdlog
   $bowtie/bowtie -p $CPU -f -n 0 -e 80 -l 18 -a -m $MAXREP --best --strata $genome $INPUT > $prj/predict.mapping.bwt
fi

echo "$perlExec $program_dir/convert_bowtie_output.pl $prj/predict.mapping.bwt > $prj/predict.mapping.arf" >> $cmdlog
$perlExec $program_dir/convert_bowtie_output.pl $prj/predict.mapping.bwt > $prj/predict.mapping.arf

echo "Parsing mapping results ..." # 2 lines
echo "$perlExec $program_dir/parse_mappings.pl $prj/predict.mapping.arf -a $MIS -b 18 -c 25 -i $MAXREP -j > $prj/predict.mapping.arf.trim" >> $cmdlog
$perlExec $program_dir/parse_mappings.pl $prj/predict.mapping.arf -a $MIS -b 18 -c 25 -i $MAXREP -j > $prj/predict.mapping.arf.trim

declare genome_seq

if [[  -e $genome.fa  ]]; then
	genome_seq=$genome.fa
elif [[ -e $genome.fas  ]]; then
	genome_seq=$genome.fas
elif [[ -e $genome.fasta  ]]; then
	genome_seq=$genome.fasta
elif [[ $genome == *".fa" || $genome == *".fas" || $genome == *".fasta" ]]; then
	genome_seq=$genome
else
	echo "Can't locate the fasta sequence for $genome";
	echo "Please provide an avaliabe $genome sequence, such as";
	echo "$genome.fa, $genome.fas, $genome.fasta, etc."
	echo "and build the Bowtie index using a command like:"
	echo "bowtie-build -f $genome.fa $genome"
	echo "exit now."
	exit 192;
fi

# 66 lines;
echo "Excise precursors candidates ..." # 7 lines
if [ $mode -eq 1 ]; then
  echo "$perlExec $program_dir/excise_precursors_iterative_final.pl $genome_seq $prj/predict.mapping.arf.trim $prj/predict.pre.fa $prj/predict.pre.coords $gear" >> $cmdlog
  $perlExec $program_dir/excise_precursors_iterative_final.pl $genome_seq $prj/predict.mapping.arf.trim $prj/predict.pre.fa $prj/predict.pre.coords $gear
else 
  echo "$perlExec $program_dir/excise_precursors_iterative_plant.pl $genome_seq $prj/predict.mapping.arf.trim $prj/predict.pre.fa $prj/predict.pre.coords $gear" >> $cmdlog
  $perlExec $program_dir/excise_precursors_iterative_plant.pl $genome_seq $prj/predict.mapping.arf.trim $prj/predict.pre.fa $prj/predict.pre.coords $gear 
fi

echo "Preparing signature  ..." # 2 lines
echo "$perlExec $program_dir/prepare_signature.pl $INPUT $prj/predict.pre.fa 1 -a $prj/$known_mature_new -o $prj/predict.signature.arf" >> $cmdlog
$perlExec $program_dir/prepare_signature.pl $INPUT $prj/predict.pre.fa 1 -a $prj/$known_mature_new -o $prj/predict.signature.arf 


echo "Fold precursors structures ..." # 2 lines
echo "$RNAfold < $prj/predict.pre.fa -noPS > $prj/predict.pre.str" >> $cmdlog
$RNAfold < $prj/predict.pre.fa -noPS > $prj/predict.pre.str
# 6 lines
if [ $rand -eq 1 ]; then
	echo "Computing randfold p-values ..."
	echo "$perlExec $program_dir/select_for_randfold.pl $prj/predict.signature.arf $prj/predict.pre.str > $prj/predict.pre.ids" >>$cmdlog
	$perlExec $program_dir/select_for_randfold.pl $prj/predict.signature.arf $prj/predict.pre.str > $prj/predict.pre.ids
	echo "$perlExec $program_dir/fastaselect.pl $prj/predict.pre.fa $prj/predict.pre.ids > $prj/predict.pre_for_rand.fa" >>$cmdlog
	$perlExec $program_dir/fastaselect.pl $prj/predict.pre.fa $prj/predict.pre.ids > $prj/predict.pre_for_rand.fa
	echo "$randfold -s $prj/predict.pre_for_rand.fa 99 > $prj/predict.pre_for_rand.rand" >>$cmdlog
	$randfold -s $prj/predict.pre_for_rand.fa 99 > $prj/predict.pre_for_rand.rand
fi

echo "Running miRDeep core algorithm ..."
declare cmd1 
declare cmd2 
declare cmd3 
if [ $mode -eq 1 ]; then
   cmd1="$perlExec $program_dir/miRDeep2_core_algorithm.pl $prj/predict.signature.arf $prj/predict.pre.str -v -50"
elif [ $mode -eq 2 ]; then
   cmd1="$perlExec $program_dir/miRDeep2_core_plant.pl $prj/predict.signature.arf $prj/predict.pre.str -v -50"
else
   cmd1="$perlExec $program_dir/miRDeep2_core_plant.pl $prj/predict.signature.arf $prj/predict.pre.str -v -50 -d"
fi

if [ $rand -eq 1 ]; then
	cmd2=`echo $cmd1`" -y $prj/predict.pre_for_rand.rand"
else
	cmd2=`echo $cmd1`
fi

if [[ $other_mature != " " ]]; then
	cmd3=`echo $cmd2`" -s $prj/$other_mature_new > $prj/predict.output.mrd"
else
	cmd3=`echo $cmd2`" > $prj/predict.output.mrd"
fi

echo $cmd3 > $prj/command_line
echo $cmd3 >>$cmdlog
sh $prj/command_line	

echo "Running permuted controls ..."
echo "$perlExec $program_dir/perform_controls.pl $prj/command_line $prj/predict.pre.str 100 $prj -a > $prj/predict.output_permuted.mrd 2>/dev/null" >>$cmdlog
$perlExec $program_dir/perform_controls.pl $prj/command_line $prj/predict.pre.str 100 $prj -a > $prj/predict.output_permuted.mrd 2>/dev/null


echo "Doing survey of accuracy ..."
echo "$perlExec $program_dir/survey.pl $prj/predict.output.mrd -a $prj/predict.output_permuted.mrd -b $prj/$known_mature_new -c $prj/predict.signature.arf -d 7 >$prj/predict.survey.csv" >>$cmdlog
$perlExec $program_dir/survey.pl $prj/predict.output.mrd -a $prj/predict.output_permuted.mrd -b $prj/$known_mature_new -c $prj/predict.signature.arf -d 7 >$prj/predict.survey.csv

echo "Make report ..."
if [[ $species != "" ]]; then
	echo "$perlExec $program_dir/make_html_simple.pl -f $prj/predict.output.mrd -p $prj/predict.pre.coords -b $score -c -e -i $prj -s $prj/predict.survey.csv -d -t $species" >>$cmdlog
	$perlExec $program_dir/make_html_simple.pl -f $prj/predict.output.mrd -p $prj/predict.pre.coords -b $score -c -e -i $prj -s $prj/predict.survey.csv -d -t $species
else
	echo "$perlExec $program_dir/make_html_simple.pl -f $prj/predict.output.mrd -p $prj/predict.pre.coords -s $prj/predict.survey.csv -i $prj -c -e -d" >>$cmdlog
	 $perlExec $program_dir/make_html_simple.pl -f $prj/predict.output.mrd -p $prj/predict.pre.coords -s $prj/predict.survey.csv -i $prj -c -e -d
fi

awk -v PROJNAME=$PROJNAME -F "\t" 'BEGIN{a=0;b=1}{if(/UCSC/){a=1}else{if(a==1){gsub("u","t",$16);gsub("u","t",$14);print ">" PROJNAME "_"b++"_"$14"_"$2"\n"$16}}}' $prj/predict.result.csv > $prj/predict.hairpin.fa
awk -v PROJNAME=$PROJNAME -F "\t" 'BEGIN{a=0;b=1}{if(/UCSC/){a=1}else{if(a==1){gsub("u","t",$14);gsub("u","t",$15);print ">"PROJNAME"_"b"_"$14"_"$2"\n"$14"\n>"PROJNAME"_"b++"_"$15"_"$2"*\n"$15}}}' $prj/predict.result.csv > $prj/predict.mature.fa
if [[ $other_mature != " " ]]; then
	echo "$perlExec $program_dir/combine_mirna_mature.pl $PROJNAME $prj/predict.mature.fa $prj/predict.hairpin.fa $prj/$other_mature_new > $prj/predict.mirna.fas" >>$cmdlog
	$perlExec $program_dir/combine_mirna_mature.pl $PROJNAME $prj/predict.mature.fa $prj/predict.hairpin.fa $prj/$other_mature_new > $prj/predict.mirna.fas
else
	echo "$perlExec $program_dir/combine_mirna_mature.pl $PROJNAME $prj/predict.mature.fa $prj/predict.hairpin.fa > $prj/predict.mirna.fas" >>$cmdlog
	$perlExec $program_dir/combine_mirna_mature.pl $PROJNAME $prj/predict.mature.fa $prj/predict.hairpin.fa > $prj/predict.mirna.fas
fi

awk -v PROJNAME=$PROJNAME -F "\t" 'BEGIN{a=0;b=1}{if(/UCSC/){a=1}else{if(a==1){gsub("u","t",$14);print PROJNAME "_"b++"_"$14"_"$2"\t"$18}}}' $prj/predict.result.csv > $prj/predict.mirna.pos

num_of_predict=`grep '>' $prj/predict.hairpin.fa | wc -l`
if [ $num_of_predict -eq 0 ] ; then
	echo "No Novel miRNA predicted, exit now."
	exit 192;
fi

echo ""
echo ""
echo "Prediction successfully done."
