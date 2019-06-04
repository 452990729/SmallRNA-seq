#!/bin/bash

basepath=$(cd `dirname $0`;pwd)

Myperl=$basepath/../BioSoft/perl
Mypython=$basepath/../BioSoft/python
raw2flN=$basepath/../BioScripts/1QC/raw2flN.pl
hq2wash=$basepath/../BioScripts/1QC/hq2wash.py
#a3='AGATCGGAAGAGCACACGTCT'
a3='TCGTATGCCGTCTTCTGCTTGT'
a5='GTTCAGAGTTCTACAGTCCGACGATC'

############clean####
$Myperl $raw2flN -i $1 -q 33 -name $2 -o $3

#############adapter####
$Mypython $hq2wash --sample $2 --qc $3/$2.stat --a3 $a3 --a5 $a5 --file $3/$2.hq.cut.fq --p 0 --outdir $3
