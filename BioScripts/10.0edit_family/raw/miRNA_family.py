# -*- coding:utf-8 -*-
#guoyang 2012.8
import re
import sys
import urllib
import urllib2
import json
import time
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(description='MiRNA family Finder')
parser.add_argument('--known',required=True,help='the known precursor')
parser.add_argument('--unknown',required=True,help="the unknown precursor")
parser.add_argument('--prefix',required=True,help="out prefix")
parser.add_argument('--family_hp',required=True,help="family hairpin to search")
parser.add_argument('--family_miFam',required=True,help="family miFam to search")
parser.add_argument('--novo_gff',required=True,help="novo miRNA hairpin sequence blast to rfam gff result")

argv=vars(parser.parse_args())
known_file=argv['known'].strip()
unknown_file=argv['unknown'].strip()
outprefix=argv['prefix'].strip()
family_hp_file=argv['family_hp'].strip()
family_miFam_file=argv['family_miFam'].strip()
novo_gff = argv['novo_gff'].strip()

#取出所有已知前体的名字
known_mir=[]
for eachLine in open(known_file):
	if eachLine.startswith('>'):
		temp=eachLine.strip().split()[0].replace('>','')
		known_mir.append(temp)

#取出所有未知前体的名字和序列
unknown_mir={}
for eachLine in open(unknown_file):
	if eachLine.startswith('>'):
		name=eachLine.replace('>','').strip()
		unknown_mir[name]=''
	else:
		unknown_mir[name]+=eachLine.strip().replace('t','u').replace('T','U').upper()
	
		
unknown_mir2={}
for line in open(novo_gff):
	if line.startswith('#'):
		continue
	else:
		temp = line.rstrip().split('\t')
		id = re.search(r'rfam-id=(\S+)',temp[-1]).group(1)
		unknown_mir2[temp[0]] = id
#构造物种简写名和全名的映射
name_map={}
for eachLine in open(family_hp_file):
	if eachLine.startswith('>'):
		temp=eachLine.strip().split()
		name_sample=re.search(r'>(.+?)-',temp[0]).group(1)
		name=' '.join(temp[2:-2])
		if name_sample not in name_map:
			name_map[name_sample]=name

#构造各个前体和家族的归属关系
family={}
for eachLine in open(family_miFam_file):
	if eachLine.startswith('ID'):
		fam_name=eachLine.strip().split()[-1]
		if fam_name not in family:
			family[fam_name]=[]
	if eachLine.startswith('MI'):
		number=eachLine.strip().split()[-1]
		family[fam_name].append(number)

species={each:OrderedDict() for each in name_map.values()}

def min_one(e):
	ord=0
	min_e=float(e[0])
	for i in range(1,len(e)):
		if float(e[i]) < min_e:
			ord=i
			min_e=float(e[i])
	return (ord,min_e)	
	



def find_family(pre):
	if pre in known_mir:
		for fa in family:
			if pre in family[fa]:
				return fa
		return ''	
	else:
		try:
			link='http://rfam.xfam.org/search/sequence'
			value={'seq':unknown_mir[pre]}
			data = urllib.urlencode(value)
			header = {'Accept': 'application/json'}
			req = urllib2.Request(link,data,header)
			response = urllib2.urlopen(req)
			obj=json.load(response)
			time.sleep(float(obj['estimatedTime'])+5)
			res_req = urllib2.Request(obj['resultURL'])
			res_req.add_header('Accept', 'application/json')
			res_response = urllib2.urlopen(res_req)
			res_obj=json.load(res_response)
			i=0
			if 'status' in res_obj:
				i+=1
				if i<20:
					print "try again..."
					time.sleep(5)
					find_family(pre)
			if res_obj['hits'] != {}:
				E=[]
				mir=[]
				for each in res_obj['hits']:
					mir.append(each)
					E.append(res_obj['hits'][each][0]['E'])
				(order,min_E)=min_one(E)
				if min_E <0.01:
					return mir[order]
				else:
					return ''
			else:
				return ''
		except:
			print pre
			return ''
			#print 'try again...'
			#find_family(pre)
		
# edit by jc on 2016/4/22	
def find_family2(pre):
	if pre in known_mir:
		for fa in family:
			if pre in family[fa]:
				return fa
		return ''
	else:
		if pre in unknown_mir2:
			if unknown_mir2[pre] in family:
				return unknown_mir2[pre]
			else:
				return ''
		else:
			return ''
#按列写入数据
mir_total=known_mir+unknown_mir.keys()
fam_list=[]
for each_pre in mir_total:
	fam=find_family2(each_pre)  #edit by jc on 2016/4/22
	print each_pre
	print fam
	for each_specie in species:
		species[each_specie][each_pre]=[]
	if (fam != '') and (family.get(fam,None)!=None):
		fam_list.append(fam)
		for each_num in family[fam]:
			name=name_map[re.match(r'(.+?)-',each_num).group(1)]
			species[name][each_pre].append(each_num)
	else:
		fam_list.append('')

#输出
output=open(outprefix+'_miRNA_family.txt','w')
output.write('\t'+'\t'.join(mir_total)+'\n')
output.write('\t'+'\t'.join(fam_list)+'\n')
for each in species:
	result=each
	for each_mir in species[each]:
		result+='\t'+' '.join(species[each][each_mir])
	output.write(result+'\n')
output.close()
	
