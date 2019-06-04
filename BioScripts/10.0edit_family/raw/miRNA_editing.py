# -*- coding:utf-8 -*-
#guoyang 2012.09.25
import sys
import os.path
import re

file=sys.argv[1].strip()
filename=os.path.basename(file)
outprefix=sys.argv[2].strip()

result=[]

def stat(sign):
	temp_result=[]
	flag3=0  #判断整体上这个mir有没有位点的改变
	k=0  #在mature上的位置
	mat_editing_count=0;
	mat_editing_count_marker={}
	start=exp.find(sign)
	end=exp.rfind(sign)
	for i in range(start, end+1):
		k+=1
		editing={}
		editing_count=0
		ori_base=pri_seq[i]
		count=0 #记录每个位点上的reads覆盖程度
		flag2=0 #记录在这个位点上是否发生了editing
		for eachsq in seq.keys():
			start1=re.search(r"([augcAUGC]+)",eachsq).start(1)
			end1=re.search(r"([augcAUGC]+)",eachsq).end(1)
			if start1 >= start and end1-1 <= end:
				count+=seq[eachsq]
				if eachsq[i] != '.' and eachsq[i].upper() != ori_base.upper():
					editing_count+=seq[eachsq]
					if not mat_editing_count_marker.has_key(eachsq):
						mat_editing_count_marker[eachsq]=seq[eachsq]
						mat_editing_count+=seq[eachsq]
					flag2=1
					flag3=1
					sub=(ori_base+'->'+eachsq[i]).upper()
					if sub in editing:
						editing[sub]+=seq[eachsq]
					else:
						editing[sub]=seq[eachsq]
		if flag2 == 1:
			temp_result+='\t\tsite'+str(k)+':\t'+str(editing_count)+'\t'+str(round(float(editing_count*100)/count,2))+'%\n'
			for each_edit in editing:
				temp_result+='\t\t\t'+each_edit+': '+str(editing[each_edit])+'\t'+str(round(float(editing[each_edit]*100)/count,2))+'%\n'
	if flag3 == 1:	
		return list(str(count)+'\t'+str(mat_editing_count)+'\t'+str(round(float(mat_editing_count*100)/count,2))+'%\n')+temp_result
	else:
		return '0\t0.00%\n'
	
for eachLine in open(file):
	if eachLine.startswith('>'):
		pre_mir=eachLine.replace('>','').strip()
		pre_mir_count=0
		pre_mir_edit=0
		mat={}
		mat_count=0
		seq={}
		flag_space=0  #记录第一个空行，这是开始统计的标志
		continue
	if eachLine.startswith('total read count'):
		temp=eachLine.strip().split()
		pre_mir_count=int(temp[-1])
		flag=1
		continue
	if eachLine.startswith('remaining read count'):
		flag=0
		continue
	if flag == 1:
		temp=eachLine.strip().split()
		mat[temp[0]]=int(temp[-1])
	if flag == 0:
		if eachLine.startswith('exp'):
			exp=eachLine.strip().split()[1]
		if eachLine.startswith('pri_seq'):
			pri_seq=eachLine.strip().split()[1]
		if re.search(r'_x\d+',eachLine):
			temp=eachLine.split()
			fold=re.search(r'_x(\d+)',temp[0]).group(1)
			if re.search(r'[AUGC]',temp[1].strip()):
				pre_mir_edit+=int(fold)
			seq[temp[1].strip()]=int(fold)
		if eachLine.strip() == '' and flag_space == 0:
			flag_space=1
			if pre_mir_count!= 0:
				result+='>pre-miRNA: '+pre_mir+'\t'+str(pre_mir_count)+'\t'+str(pre_mir_edit)+'\t'+str(round(float(pre_mir_edit*100)/pre_mir_count,2))+'%\n'
				for each in mat:
					if mat[each] != 0:
						mat_count=mat[each]
						result+='\tmature miRNA: '+each+'\t'
						if re.search(r'-5p$',each):
							result+=stat('5')
						elif re.search(r'-3p$',each):
							result+=stat('3')
						elif re.search(r'\*$',each):
							result+=stat('*')
						else:
							result+=stat('M')
			else:
				result+='>pre-miRNA: '+pre_mir+'\t0\t0\t0.00%\n'
		
open(outprefix+'_editing_stats.txt','w').writelines(result)		

