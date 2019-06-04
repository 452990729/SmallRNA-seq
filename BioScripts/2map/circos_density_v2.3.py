import os
import sys
import HTSeq
import random
import os.path
import argparse
import ConfigParser
from Bio import SeqIO
import rpy2.robjects as robjects
#guoyang 20130723

parser = argparse.ArgumentParser(description='coverage density distribution -- circos guoyang@novogene.cn')
parser.add_argument('--bam',required=True,help="the bam file from mapping")
parser.add_argument('--fa',required=True,help="genome fasta file")
parser.add_argument('--n',type=int,default=10,help="n chrom or contig will be drew default:10")
parser.add_argument('--iv',type=int,default=1000,help="the interval for coverage default:1000")
parser.add_argument('--r',type=int,default=10000,help="r reads for reads distribution track default:10000")
parser.add_argument('--outputfile',default=None,help="Change the output filename.")
parser.add_argument('--outputdir',default='./',help="Output directory")
args=parser.parse_args()
bam=args.bam
assert os.path.isfile(bam)
fa=args.fa
assert os.path.isfile(fa)
n=args.n
iv=args.iv
r=args.r
outputfile = args.outputfile
outdir  = args.outputdir

mean=robjects.r['mean']
sd=robjects.r['sd']

script_path = os.path.dirname(os.path.abspath(__file__))
config = ConfigParser.ConfigParser()
config.read('{}/../../BioModule/Settings/config.ini'.format(script_path))
circos_path = config.get('software','circos')
perlExec = config.get('srnaenv','perl_v5182')

#genome file length
genome_lenth={}
chrom_count=0
for seq_record in SeqIO.parse(fa,'fasta'):
	chrom_count+=1
	if seq_record.id not in genome_lenth:
		genome_lenth[seq_record.id]=len(seq_record.seq)
	else:
		print "please attention,genome seq_record.id duplicate"
print genome_lenth
#n=min(n,chrom_count)

#when genome include chr,must select chr level

#sort genome_lenth by lenth from big to small 
#sorted_genome_lenth=reversed(sorted(genome_lenth.values()))
sorted_genome_lenth=sorted(genome_lenth.items(),key=lambda genome_lenth:genome_lenth[1],reverse=True)
lenth_ten={}

if "chr1" in genome_lenth.keys():
	print "genome is chr level,select chr level for circos"
	num=0
	for i in range(1,11):
		chr_name="chr"+str(i)
		num+=1
		if chr_name in genome_lenth.keys():
			lenth_ten[chr_name]=genome_lenth[chr_name]
		else:
			break
elif "1" in genome_lenth.keys():
	print "genome is chr level,select chr level for circos"
        num=0
	for i in range(1,11):
		chr_name=str(i)
		num+=1
		if chr_name in genome_lenth.keys():
			lenth_ten[chr_name]=genome_lenth[chr_name]
		else:
			break
else:
	print "genome is not chr level,select ten most longest ten scafoold or contig  for circos"
	if len(genome_lenth.keys()) >9:
		for i in range(1,11):
			lenth_ten[sorted_genome_lenth[i][0]]=sorted_genome_lenth[i][1]
		num=11
	else:
		for i in range(1,len(genome_lenth.keys()+1)):
			lenth_ten[sorted_genome_lenth[i][0]]=sorted_genome_lenth[i][1]
		num=len(genome_lenth.keys())+1
chr_num=num-1
n=min(chr_num,chrom_count)
total_len=sum(lenth_ten.values())
if total_len > 3116677:
	units=int(round(float(total_len)/3116677))*10000
else:
	units=10000

genome=[]
flag=0
for each in lenth_ten:
	if flag == 0:
		genome.append('chr - '+each+' '+each+' 0 '+str(genome_lenth[each]-1)+' vlgrey\n')
		flag=1
		continue
	if flag == 1:
		genome.append('chr - '+each+' '+each+' 0 '+str(genome_lenth[each]-1)+' grey\n')
		flag=0
open(os.path.join(outdir,'genome.txt'),'w+').writelines(genome)	


#bam file
bam_file=HTSeq.BAM_Reader(bam)
coverage = HTSeq.GenomicArray( "auto", stranded=True, typecode="i" )
reads=[]
read_count=0
for almnt in bam_file:
	if almnt.aligned and almnt.iv.chrom in lenth_ten:
		coverage[almnt.iv]+=1
		read_count+=1
		if almnt.iv.strand == '+':
			reads.append(almnt.iv.chrom+'\t'+str(almnt.iv.start)+'\t'+str(almnt.iv.end-1)+'\tstroke_color=vdred\n')
		if almnt.iv.strand == '-':
			reads.append(almnt.iv.chrom+'\t'+str(almnt.iv.start)+'\t'+str(almnt.iv.end-1)+'\tstroke_color=vdblue\n')


#read tile distribution
r=min(r,read_count)
read_sample=random.sample(reads,r)
open(os.path.join(outdir,'reads.txt'),'w+').writelines(read_sample)


#line density
line_plus=[]
line_minus=[]
coverage_plus_total=[]
coverage_minus_total=[]
for each in lenth_ten:
	flag=1
	start=-iv
	while True:
		start+=iv
		if start+iv <= lenth_ten[each]:
			window_plus=HTSeq.GenomicInterval(each,start,start+iv,'+')
			coverage_plus=list(coverage[window_plus])
			coverage_plus_mean=mean(robjects.IntVector(coverage_plus))[0]
			coverage_plus_total.append(coverage_plus_mean)
			line_plus.append(each+'\t'+str(start)+'\t'+str(start+iv-1)+'\t'+str(coverage_plus_mean)+'\n')
			window_minus=HTSeq.GenomicInterval(each,start,start+iv,'-')
			coverage_minus=list(coverage[window_minus])
			coverage_minus_mean=mean(robjects.IntVector(coverage_minus))[0]
			coverage_minus_total.append(coverage_minus_mean)
			line_minus.append(each+'\t'+str(start)+'\t'+str(start+iv-1)+'\t'+str(coverage_minus_mean)+'\n')
		else:
			window_plus=HTSeq.GenomicInterval(each,start,lenth_ten[each],'+')
			coverage_plus=list(coverage[window_plus])
			coverage_plus_mean=mean(robjects.IntVector(coverage_plus))[0]
			coverage_plus_total.append(coverage_plus_mean)
			line_plus.append(each+'\t'+str(start)+'\t'+str(lenth_ten[each]-1)+'\t'+str(coverage_plus_mean)+'\n')
			window_minus=HTSeq.GenomicInterval(each,start,lenth_ten[each],'-')
			coverage_minus=list(coverage[window_minus])
			coverage_minus_mean=mean(robjects.IntVector(coverage_minus))[0]
			coverage_minus_total.append(coverage_minus_mean)
			line_minus.append(each+'\t'+str(start)+'\t'+str(lenth_ten[each]-1)+'\t'+str(coverage_minus_mean)+'\n')
			flag=0
		if flag == 0:
			break

#filter
plus_cutoff=mean(robjects.FloatVector(coverage_plus_total))[0]+3*sd(robjects.FloatVector(coverage_plus_total))[0]
minus_cutoff=mean(robjects.FloatVector(coverage_minus_total))[0]+3*sd(robjects.FloatVector(coverage_minus_total))[0]		
line_plus_filtered=[]
line_minus_filtered=[]
for eachLine in line_plus:
	if float(eachLine.split()[-1].strip()) < plus_cutoff:
		line_plus_filtered.append(eachLine)
for eachLine in line_minus:
	if float(eachLine.split()[-1].strip()) < minus_cutoff:
		line_minus_filtered.append(eachLine)	

open(os.path.join(outdir,'line_plus.txt'),'w').writelines(line_plus_filtered)
open(os.path.join(outdir,'line_minus.txt'),'w').writelines(line_minus_filtered)


ideogram_position='''
radius           = 0.775r
thickness        = 30p
fill             = yes
fill_color       = black
stroke_thickness = 2
stroke_color     = black
'''
open(os.path.join(outdir,'ideogram.position.conf'),'w+').write(ideogram_position)

ticks='''

show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = dims(ideogram,radius_outer)
orientation      = out
label_multiplier = 1e-4
color            = black

thickness        = 2p
label_offset     = 5p
format           = %d

<tick>
spacing        = 1u
show_label     = no
size           = 10p
</tick>

<tick>
spacing        = 5u
show_label     = yes
label_size     = 15p
size           = 15p
</tick>

<tick>
spacing        = 10u
show_label     = yes
label_size     = 20p
size           = 20p
</tick>

</ticks>

'''
open(os.path.join(outdir,'ticks.conf'),'w+').write(ticks)


ideogram='''

<ideogram>

<spacing>
default = 0.01r
#break   = 0.5r
</spacing>

<<include ideogram.position.conf>>
<<include ideogram.label.conf>>


radius*       = 0.75r

</ideogram>

'''
open(os.path.join(outdir,'ideogram.conf'),'w+').write(ideogram)


ideogram_label='''
show_label       = yes
label_font       = default
label_radius     = dims(ideogram,radius) + 0.1r
label_size       = 24
label_parallel   = yes
label_case       = lower

'''
open(os.path.join(outdir,'ideogram.label.conf'),'w+').write(ideogram_label)

circos='''

<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>



karyotype   = genome.txt

<image>
<<include etc/image.conf>>
</image>

chromosomes_units = %s

chromosomes_display_default = yes

<plots>

<plot>
type      = tile
file      = reads.txt
r1        = 0.95r
r0        = 0.75r
layers    = 4
margin    = 0.2u
thickness = 4
padding   = 2
layers_overflow  = grow
orientation      = in
stroke_thickness = 1p
#stroke_color     = grey
#color            = orange

<backgrounds>
<background>
color     = vvlgrey
#y0        = 0.006
</background>

</backgrounds>

</plot>


<plot>
type        = histogram
file        = line_plus.txt
r1          = 0.70r
r0          = 0.55r
#min         = 0
#max         = 1
extend_bin  = yes
fill_color  = vdorange
color       = vdorange
thickness   = 0
orientation = out

<axes>
<axis>
color     = lgrey_a2
thickness = 1
spacing   = 0.1r
</axis>
</axes>

</plot>

<plot>
type        = histogram
file        = line_minus.txt
r1          = 0.55r
r0          = 0.40r
#max         = 1
#min         = 0
extend_bin  = yes
fill_color  = vdgreen
color       = vdgreen
thickness   = 0
orientation = in

<axes>
<axis>
color     = lgrey_a2
thickness = 1
spacing   = 0.1r
</axis>
</axes>
</plot>


</plots>

<<include etc/housekeeping.conf>>

''' % (units)
open(os.path.join(outdir,'circos.conf'),'w+').write(circos)

if outputfile:
	os.system('%s %s/bin/circos -conf %s/circos.conf --outputdir %s --outputfile %s' % (perlExec,circos_path,outdir,outdir,outputfile))
else:
	os.system('%s %s/bin/circos -conf %s/circos.conf --outputdir %s ' % (perlExec,circos_path,outdir,outdir))
