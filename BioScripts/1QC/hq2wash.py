# -*- coding:utf-8 -*-
import sys
import argparse
import re
import os.path
import Levenshtein
import rpy2.robjects as robjects

parser = argparse.ArgumentParser(description='adapter(3 and 5) removal')
parser.add_argument('--sample',required=True,help='Sample Name')
parser.add_argument('--qc',required=True,help='the stat report after QC, including total, N%%>10%%, low qual')
parser.add_argument('--a3',help='the adapter 3 sequence, default: CTGTAGGC',default='CTGTAGGC')
parser.add_argument('--a5',help="the adapter 5 sequence, default: CGACGATC",default='CGACGATC')
parser.add_argument('--p',type=int,help="cut base number from 5', default:3",default=3)
parser.add_argument('--min3',type=int,help="the minimum length for adapter3, default: 8",default=8)
parser.add_argument('--min5',type=int,help="the minimum length for adapter5, default: 8",default=8)
parser.add_argument('--file',required=True,help="the read fastq file")
parser.add_argument('--min_len',type=int,help="the minimum length of reads after filtered, default: 16",default=16)
parser.add_argument('--max_len',type=int,help="the maximum length of reads after filtered, default: 35",default=35)
parser.add_argument('--outdir', metavar = 'path', help = 'The output files path',required=True)
argv=vars(parser.parse_args())


outdir = argv['outdir'].strip()
sample=argv['sample'].strip()
adapter3=argv['a3'].strip()
adapter5=argv['a5'].strip()
part=argv['p']
the_minimum_length_of_adapter3=argv['min3']
the_minimum_length_of_adapter5=argv['min5']
file=argv['file'].strip()
filename=os.path.basename(file)
min_len=argv['min_len']
max_len=argv['max_len']

if len(adapter3) <  the_minimum_length_of_adapter3:
	sys.exit('the 3 end adapter sequence given is too short... given:%s' % the_minimum_length_of_adapter3)
if len(adapter5) <  the_minimum_length_of_adapter5:
	sys.exit('the 5 end adapter sequence given is too short... given:%s' % the_minimum_length_of_adapter5)


def trim5(seq):
	for base in range(the_minimum_length_of_adapter5,len(seq)):#取不到的那个值
		if len(seq[:base]) < len(adapter5):
			seq_part=seq[:base]
			ada_part=adapter5[-base:]
			if Levenshtein.hamming(seq_part,ada_part) <= 1:
				return True
		else:
			seq_part=seq[base-len(adapter5):base]
			ada_part=adapter5
			if Levenshtein.hamming(seq_part,ada_part) <= 1:
				return True
	return False	

def trim3(seq):
	read_trimed=''
	for base in range(len(seq)-the_minimum_length_of_adapter3,-1,-1):
		if len(seq[base:]) < len(adapter3):
			seq_part=seq[base:]
			ada_part=adapter3[:len(seq_part)]
			if Levenshtein.hamming(seq_part,ada_part) <= 2:
				read_trimed=seq[:base]
		else:
			seq_part=seq[base:base+len(adapter3)]
			ada_part=adapter3
			if Levenshtein.hamming(seq_part,ada_part) <=2:
				read_trimed=seq[:base]
	return read_trimed

A=re.compile(r'A{7,}')
T=re.compile(r'T{7,}')
G=re.compile(r'G{7,}')
C=re.compile(r'C{7,}')

def test_poly(seq):
	return bool(A.search(seq)) or bool(T.search(seq)) or bool(G.search(seq)) or bool(C.search(seq))
			
def test_length(seq):
	if len(seq)>=min_len and len(seq)<=max_len:
		return True
	else:
		return False

temp_all=[]	
for eachLine in open(argv['qc']):
	if eachLine.strip() == '':
		break
	temp=eachLine.strip().split('\t')
	temp_all.append(temp[1].strip())

total=int(temp_all[1])
total_base=int(temp_all[2])
read_len=total_base/total
polyN=int(temp_all[10])
polyN_fbase=polyN*read_len
low_qual=int(temp_all[11])
low_qual_fbase=low_qual*read_len

sequence={}
with5_count=0
with5_fbase=0
insert_illegal_count=0	#without 3 adapter or len(insert sequences) < 10
insert_illegal_fbase=0
trim3_base=0
poly_count=0
poly_fbase=0
clean=0
clean_base=0
length_count=0
length_fbase=0
remain=0
remain_base=0

counter=0
clean_total_fasta=open(os.path.join(outdir,sample+'_clean_total.fa'),'w+')	#clean_total_fasta=[]
remain_total_fasta=open(os.path.join(outdir,sample+'_remain_total.fa'),'w+')	#remain_total_fasta=[]
reads_file=open(file,'r').readlines()
for eachLine in reads_file:
	if eachLine.strip() != '':
		counter+=1
		if counter == 4:
			counter=0
		if counter == 1:
			read_name='>'+eachLine
		if counter == 2:
			temp=eachLine.strip()
			if trim5(temp) == True: #有5端adapter
				with5_count+=1
				with5_fbase+=len(temp)+part
				continue
			temp_trim3=trim3(temp)
			if len(temp_trim3) < 10: #没有3端adapter
				insert_illegal_count+=1
				insert_illegal_fbase+=len(temp)+part
				continue
			trim3_base+=len(temp)-len(temp_trim3)+part
			if test_poly(temp_trim3):#有polyATGC
				poly_count+=1
				poly_fbase+=len(temp_trim3)
				continue
			clean+=1
			clean_base+=len(temp_trim3)
			clean_total_fasta.write(read_name)	#clean_total_fasta.append(read_name)
			clean_total_fasta.write(temp_trim3+'\n')	#clean_total_fasta.append(temp_trim3+'\n')
			if not test_length(temp_trim3):
				length_count+=1
				length_fbase+=len(temp_trim3)
				continue
			remain+=1
			remain_base+=len(temp_trim3)
			remain_total_fasta.write(read_name)	#remain_total_fasta.append(read_name)
			remain_total_fasta.write(temp_trim3+'\n')	#remain_total_fasta.append(temp_trim3+'\n')
			if temp_trim3 not in sequence:
				sequence[temp_trim3]=1
			else:
				sequence[temp_trim3]+=1

clean_total_fasta.close()
remain_total_fasta.close()

remain_uniq_fasta=open(os.path.join(outdir,sample+'_remain_uniq.fa'),'w+')		
#remain_uniq_fasta=[]
remain_uniq=0
remain_uniq_base=0
for each in sequence:
	remain_uniq+=1
	remain_uniq_base+=len(each)
	remain_uniq_fasta.write('>'+each+'('+str(sequence[each])+')'+'\n')
	remain_uniq_fasta.write(each+'\n')
	#remain_uniq_fasta.append('>'+each+'('+str(sequence[each])+')'+'\n')
	#remain_uniq_fasta.append(each+'\n')
remain_uniq_fasta.close()
#open(sample+'_clean_total.fa','w').writelines(clean_total_fasta)
#open(sample+'_remain_total.fa','w').writelines(remain_total_fasta)
#open(sample+'_remain_uniq.fa','w').writelines(remain_uniq_fasta)

overview=open(os.path.join(outdir,sample+'_clean_process_overview.txt'),'w+')
overview.write('total reads:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (total,100,'%',total_base,100,'%'))
overview.write('filter reads related to N:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (polyN,float(polyN*100)/total,'%',polyN_fbase,float(polyN_fbase*100)/total_base,'%'))
overview.write('filter reads related to low qual:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (low_qual,float(low_qual*100)/total,'%',low_qual_fbase,float(low_qual_fbase*100)/total_base,'%'))
overview.write('filter reads with 5 end adapter:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (with5_count,float(with5_count*100)/total,'%',with5_fbase,float(with5_fbase*100)/total_base,'%'))
overview.write('filter reads 3_adapter_null or insert_null:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (insert_illegal_count,float(insert_illegal_count*100)/total,'%',insert_illegal_fbase,float(insert_illegal_fbase*100)/total_base,'%'))
overview.write('trim reads 3 end adapter:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (0,0,'%',trim3_base,float(trim3_base*100)/total_base,'%'))
overview.write('filter reads with polyATGC:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (poly_count,float(poly_count*100)/total,'%',poly_fbase,float(poly_fbase*100)/total_base,'%'))
overview.write('clean reads:\t%s\t(%s %s)\t%s\t(%s %s)\n\n' % (clean,float(clean*100)/total,'%',clean_base,float(clean_base*100)/total_base,'%'))

overview.write('filter reads with illegal length:\t%s\t(%s %s)\t%s\t(%s %s)\n' % (length_count,float(length_count*100)/total,'%',length_fbase,float(length_fbase*100)/total_base,'%'))
overview.write('remined reads after length filter:\t%s\t(%s %s)\t%s\t(%s %s)\n\n' % (remain,float(remain*100)/total,'%',remain_base,float(remain_base*100)/total_base,'%'))

overview.write('Parameters:\n')
overview.write('the read fastq file: %s\n' % filename)
overview.write('the adapter 3 sequence: %s\n' % adapter3)
overview.write('the minimum length for adapter3: %s\n' % the_minimum_length_of_adapter3)
overview.write('the adapter 5 sequence: %s\n' % adapter5)
overview.write('the minimum length for adapter5: %s\n' % the_minimum_length_of_adapter5)
overview.write('the minimum length of reads: %s\n' % min_len)
overview.write('the maximum length of reads: %s\n' % max_len)
overview.close()

overview=open(os.path.join(outdir,sample+'_remain_total_uniq.txt'),'w+')
overview.write('total reads:\t%s\n' % remain)
overview.write('total bases (bp):\t%s\n' % remain_base)
overview.write('uniq reads:\t%s\n' % remain_uniq)
overview.write('uniq bases (bp):\t%s\n' % remain_uniq_base)
overview.close()

robjects.r('pieval=c(%s)' % ','.join(map(str,[polyN,low_qual,with5_count,insert_illegal_count,poly_count,clean])))
robjects.r('pielabels=c(%s)' % ("'"+"','".join(['N%% > 10%%(%s, %s%%)' % (polyN,round(float(polyN*100)/total,2)),'low quality(%s, %s%%)' % (low_qual,round(float(low_qual*100)/total,2)),'with 5 end adapter(%s, %s%%)' % (with5_count,round(float(with5_count*100)/total,2)),'3_adapter_null or insert_null(%s, %s%%)' % (insert_illegal_count,round(float(insert_illegal_count*100)/total,2)),'with polyA/T/G/C(%s, %s%%)' % (poly_count,round(float(poly_count*100)/total,2)),'clean reads(%s, %s%%)' % (clean,round(float(clean*100)/total,2))])+"'"))
robjects.r('sample=c(%s)' % ("'"+sample+"'"))
robjects.r('outdir=c(%s)' % ("'"+outdir+"'"))
robjects.r('''
library('ggplot2')
pieval[pieval==0]<-NA
dat=data.frame(pieval,pielabels)
p<-ggplot(na.omit(dat),aes(x='',y=pieval,fill=pielabels))+geom_bar(width = 1,stat="identity")+scale_y_continuous(breaks=c())+coord_polar(theta = "y")+theme(legend.title=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank())+ggtitle(paste("Classification of raw data\n","(", sample,")",sep=""))+scale_fill_brewer(palette="Set1")
ggsave(filename=paste(outdir,'/',sample,"_clean_process_overview_pie.png",sep=""),type='cairo',dpi=300)
dev.off()
ggsave(filename=paste(outdir,'/',sample,"_clean_process_overview_pie.pdf",sep=""))
dev.off()
''')

if os.path.exists('Rplots.pdf'):
    os.remove('Rplots.pdf')
