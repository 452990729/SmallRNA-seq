#!/usr/bin/env python
# -*- coding=utf-8 -*-

__author__  = 'wangyunkai'
__date__    = '2018-3-30'

import re
import os
import glob
import linecache
from itertools import chain
from collections import defaultdict

class CheckTools(object):

    def __init__(self,args):
        self.args          =  args
        self.projpath      =  args.get('projpath') or os.getcwd()
        self.runshell      =  os.path.join(self.projpath,'run.sh')
        self._samples      =  args.get('sample').strip().split(',')       
        self._groupname    =  args.get('groupname')
        self._group        =  args.get('group')
        self._compare      =  args.get('compare')
        self._venn_cluster =  args.get('venn_cluster')
        self.org           =  args.get('org')
        self.min_len       =  18
        self._max_len      =  30
        if args.get('mapfile'):
            self.mapfile   =  os.path.join(self.projpath,args.get('mapfile'))
        self.mdspe         =  args.get('mdspe')
        self.mdspedb       =  args.get('mdspedb')
        self.refer         =  args.get('refer')
        self.rRNA          =  args.get('rRNA')
        self.tRNA          =  args.get('tRNA')        
        self.snRNA         =  args.get('snRNA')
        self.snoRNA        =  args.get('snoRNA')
        self.repdir        =  args.get('repdir')
        self.gtf           =  args.get('gtf')
        self.exon          =  args.get('exon')
        self.intron        =  args.get('intron')
        self.gene          =  args.get('gene')
        self.geneAnn       =  args.get('geneAnn')
        self.go            =  args.get('go')
        self.utr3          =  args.get('utr3')
 

    _TUE  = '\033[36m'
    _FLE  = '\033[31m'
    _ENDC = '\033[0m'
    _WARNING = '\033[33m'

    def checkrunshell(self):
        '''
        check runshell config file
        '''
        filelist = [ each.strip('\n\ ').split()[-1] for each in \
                   linecache.getlines(os.path.join(self.projpath,self.runshell)) if '/' in each and 'sched' not in each]
        for eachfile in filelist:
            if os.path.exists(eachfile):
                print "%s %s%s%s" %(eachfile,self._TUE,os.path.exists(eachfile),self._ENDC)
            else:
                print "%s %s%s%s" %(eachfile,self._FLE,os.path.exists(eachfile),self._ENDC)


    @property
    def max_len(self):    
        if 'refplant' in self.org:
            return self._max_len
        else:
            self._max_len = 35
            return self._max_len
 
    @property
    def samples(self):
        '''
            Naming Notations
        '''
        message = '''can only use word/digit/underscore in sample, start with word. 
        Length is 8, no windows reserved words like CON PRN...'''
        PatChar = r'^(CON|PRN|AUX|CLOCK\$|NUL|COM[1-9]|LPT1|[\W]+)$'
        PatNum = r'^\d+'
        MAXLEN = 50
        for sampid in self._samples:
            if re.search(PatChar,sampid) or re.search(PatNum,sampid) or len(sampid) > MAXLEN:
                print self._FLE,'\tinvalid sample name==>{}'.format(sampid),self._ENDC
                print self._WARNING,message,self._ENDC
                exit()

        '''
        样本名称不能通过大小写来区分
        '''
        if len(self._samples) != len(set([i.upper() for i in self._samples])):
            print self._FLE, "\t样本名称不能通过大小写来区分\t",self._ENDC
            exit()
        return self._samples 
 

    @property
    def group(self):

        if self._group == None: 
            groups        = self.samples          
        else:
            groups = [
                      each.strip().split(':') for each in self._group.strip(' :').split(',') \
                      if each.strip() != ''
                     ]
        self._group=','.join([':'.join(each) for each in groups])
              
        if len(groups) == 1:
            sample_list = groups
        else:      
            sample_list = list(chain.from_iterable(groups))
        assert set(sample_list).issubset(self.samples)        
        return self._group


    @property
    def groupname(self): 
        '''
        分组名称不能通过大小写来区分
        '''
        if self._groupname:
            if len(self._groupname.strip().split(',')) != len(set([i.upper() for i in self._groupname.strip().split(',')])):
                print self._WARNING, "\t组合名称不能通过大小写来区分\t",self._ENDC
                exit()
        
        global groupnames
          
        groups = self.group.strip(' :').split(',')
        if self._groupname == None:
            if ':' in groups:                                 #含有生物学重复
                groupnames = [each.split(':') for each in groups]
            else:
                groupnames = ['group'+str(k+1) for k in range(len(groups))]     
        else:
            groupnames = [ each.strip() for each in self._groupname.strip().split(',') if each.strip() != '']
 
        assert len(groupnames) == len(groups)
        self._groupname = ','.join(groupnames)
        return self._groupname


    @property
    def compare(self):

        groupname = self.groupname
        if self._compare != None:
            groupname_dict = dict()
            compares_temp  = list()
            if os.path.isfile(self._compare):       #compare 文件格式是 组合1 groupname1  groupname2 中间以空格或者\t分隔
                for g_index,g in enumerate(groupnames):
                    groupname_dict[g] = g_index + 1
                with open(self._compare) as f:
                    for line in f:
                        temp = line.split()
                        compare_temp = '{}:{}'.format(groupname_dict[temp[1]],groupname_dict[temp[2]])
                        compares_temp.append(compare_temp)
                compare = ','.join(compares_temp)
            else:
                compare = self._compare
        else:
                compare='1:1'
        compares=[each.strip().split(':') for each in compare.strip().split(',') if each.strip() != '']

        group_count = list()
        total_group_name = list()
        for each_compare in compares:
            assert len(each_compare) == 2
            single_group_name = list()
            for each_group in each_compare:
                assert each_group.isdigit()
                group_count.append(int(each_group))
                single_group_name.append(groupnames[int(each_group)-1])
            total_group_name.append(':'.join(single_group_name))
        assert max(group_count) <= len(groupnames)
        assert min(group_count) > 0   
        
        compare_name = ','.join(total_group_name)
        compare = ','.join([':'.join(each) for each in compares])
        self._compare = compare
        return self._compare,compare_name


    @property
    def venn_cluster(self):

        groupnames = [ each.strip() for each in self.groupname.split(',') ]
        if self._venn_cluster != None:
            if os.path.isfile(self._venn_cluster):
                compare_temp= self.compare[0].split(',')
                ven_clust = []
                with open(self._venn_cluster) as f:           
                    for line in f:
                        element = []
                        temp = re.findall(r'\d+',line)[1:]
                        transfm = [int(i)-1 for i in temp]
                        for i in transfm:
                            element.append(compare_temp[i])
                        element = '_'.join(element)
                        ven_clust.append(element)
                venn_cluster = ','.join(ven_clust)
            else:
                venn_cluster=self._venn_cluster.strip()
            com_pairs=self.compare[0].split(',')
            venn_clusters=[ each.split('_') for each in venn_cluster.split(',')]
            temp1=[]
            for each1 in venn_clusters:
                temp2=[]
                for each2 in each1:
                    assert each2 in com_pairs
                    temp3=each2.split(':')
                    assert len(temp3) == 2
                    temp2.append(groupnames[int(temp3[0])-1]+':'+groupnames[int(temp3[1])-1])
                temp1.append('_'.join(temp2))   
            venn_cluster_name=','.join(temp1)   #venn_cluster_name  
        else:
            venn_cluster = None
            venn_cluster_name = None
        self._venn_cluster = venn_cluster
        return self._venn_cluster,venn_cluster_name


    def checkSErawdata(self):
        '''
        检查rawdata是否需要合并
        '''
        rawfq_dict = defaultdict(set)
        rawadapter_dict = defaultdict(set)
        with open(self.mapfile) as file:
            for line in file:
                array = line.strip().split('\t')
                raw_path = array[0]
                library  = array[1]
                sampid   = array[2]
                for fq in glob.glob(os.path.join(raw_path,library,'{}*.fq.gz'.format(library))):
                    rawfq_dict[sampid].add(fq)
                for adapter in glob.glob(os.path.join(raw_path,library,'{}*.adapter.list.gz'.format(library))):
                    rawadapter_dict[sampid].add(adapter)
                assert len(rawfq_dict[sampid]) == len(rawadapter_dict[sampid]), \
                      (self._WARNING,'{sampid}样本的rawdata或者adapter有一个不存在'.format(sampid = sampid),self._ENDC)
        return rawfq_dict,rawadapter_dict


    def getRefConfig(self):
        '''
        获取基因组配置文件
        '''
        dict_refconfig = defaultdict(dict)
        if self.mdspe:
            mdspe = self.mdspe.strip()
            genome_dir='{}'.format(self.mdspedb)+'/'+mdspe[0:3]+'/'+mdspe[4:]+'/'+mdspe
            dict_refconfig['refer']   = self.checkfile(genome_dir+'.fa') 
            dict_refconfig['rRNA']    = self.checkfile(genome_dir+'_rRNA.fa')
            dict_refconfig['tRNA']    = self.checkfile(genome_dir+'_tRNA.fa')
            dict_refconfig['snRNA']   = self.checkfile(genome_dir+'_snRNA.fa')
            dict_refconfig['snoRNA']  = self.checkfile(genome_dir+'_snoRNA.fa')
            dict_refconfig['repdir']  = self.checkfile(genome_dir+'_repeat_predict')
            dict_refconfig['gtf']     = self.checkfile(genome_dir+'.gtf')
            dict_refconfig['geneAnn'] = self.checkfile(genome_dir+'_genenamefile')
            dict_refconfig['go']      = self.checkfile(genome_dir+'.go')
            if self.org.lower() == 'refplant':
                dict_refconfig['exon']   = self.checkfile(genome_dir+'_exon.fa')
                dict_refconfig['intron'] = self.checkfile(genome_dir+'_intron.fa')
                dict_refconfig['gene']   = self.checkfile(genome_dir+'_transcript.fa')
            elif self.org.lower() == 'refanimal':
                dict_refconfig['utr3']   = self.checkfile(genome_dir+'_transcript_3_utr.fa')
                if 'hsa' in self.mdspe.lower():
                    genome_dir='{}'.format(self.mdspedb)+'/'+mdspe[0:3]+'/'+mdspe[4:]
                    dict_refconfig['exon']   = self.checkfile(genome_dir+'/exon/'+ mdspe+'_exon.fa')
                    dict_refconfig['intron'] = self.checkfile(genome_dir+'/intron/'+mdspe+'_intron.fa')
                else:
                    dict_refconfig['exon']   = self.checkfile(genome_dir+'_exon.fa')
                    dict_refconfig['intron'] = self.checkfile(genome_dir+'_intron.fa')
        else:
            dict_refconfig['refer']   = self.refer   
            dict_refconfig['rRNA']    = self.rRNA
            dict_refconfig['tRNA']    = self.tRNA
            dict_refconfig['snRNA']   = self.snRNA
            dict_refconfig['snoRNA']  = self.snoRNA
            dict_refconfig['repdir']  = self.repdir
            dict_refconfig['gtf']     = self.gtf
            dict_refconfig['exon']    = self.exon
            dict_refconfig['intron']  = self.intron
            dict_refconfig['gene']    = self.gene
            dict_refconfig['geneAnn'] = self.geneAnn
            dict_refconfig['go'] = self.go
            dict_refconfig['utr3'] = self.utr3

        return dict_refconfig


    def checkfile(self,file):
        if os.path.exists(file):
            return file
        else:
            print self._FLE,'file not exists :%s ' % file ,self._ENDC
            exit()
