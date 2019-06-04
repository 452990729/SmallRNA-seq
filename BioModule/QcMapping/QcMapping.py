#!/usr/bin/env python
# _*_ coding:utf-8 _*_


import os
import glob
from Settings import BioConfig
from itertools import combinations
from collections import defaultdict
from Stools import CheckTools, AtomJob, Job, write_shell


class QcMapping(object):
    
    def __init__(self,args,job):
        self.args         = args
        self.projpath     = args.get('projpath') or os.getcwd()
        self.mapfile      = os.path.join(self.projpath,args.get('mapfile'))
        self.org          = args.get('org')
        self.sched        = args.get('sched')
        self.project      = args.get('project')
        self.chrNum       = args.get('chrNum')
        self.refer        = CheckTools(self.args).getRefConfig().get('refer')
        self.sample_list  = CheckTools(self.args).samples
        self.sample       = args.get('sample')
        self.common       = args.get('common')
        self.ad3          = args.get('ad3')
        self.ad5          = args.get('ad5')
        self.min_len      = CheckTools(self.args).min_len
        self.max_len      = CheckTools(self.args).max_len
        self.analycode    = args.get('analy_code') if args.get('analy_code') else BioConfig.ORG_ANALYSIS_CODE[self.org]
        self.job          = job
        self.BioScripts   = BioConfig.config.get('srnaenv','BioScripts')
        self.perlExec     = BioConfig.config.get('srnaenv','perl_v5240')
        self.R_v2153      = BioConfig.config.get('srnaenv','R_v2153')
        self.python276    = BioConfig.config.get('srnaenv','python_v276')
        self.LIBRARY_PATH = BioConfig.config.get('srnaenv','LIBRARY_PATH')
        self.R_v303       = BioConfig.config.get('srnaenv','R_v303')
        self.bowtie       = BioConfig.config.get('software','bowtie1')
        self.bowtie2      = BioConfig.config.get('software','bowtie2')
        self.circos       = BioConfig.config.get('software','circos')
        self.samtools     = BioConfig.config.get('software','samtools')
        self.abbr         = args.get('abbr')
        if args.get('common') == 'y':
            self.qc_code  = 1.2
        else:
            if '1.1' in self.analycode:
                self.qc_code = 1.1
            else:
                self.qc_code = 0

        if 'noref' in self.org:
            self.map_code = 2.1      #不画circos图
        else:
            self.map_code = 2.2 
 

    def run(self):
        
        CheckTools(self.args).checkrunshell()
        for sampid in self.sample_list:
            self.draw_length_distribution(sampid)
            self._category(sampid)
            self._MapStat(sampid)
            if self.map_code == 2.2:
                self._circos(sampid)
        #self._qc_report()
        if self.qc_code == 1.2:
            for group in [c for c in  combinations(self.sample_list, 2)]:
                self._sharespecific(group)
                
        
    def draw_length_distribution(self,sampid):

        dict_clean_sample = dict()
        with open(self.mapfile) as mf:
            for line in mf:
                if line:
                    array = line.strip().split('\t')
                    dict_clean_sample[array[1]] = array[0]
                
        code = 'date +"%D %T -> Start Stat sRNA distribution " && \\\n'

        if dict_clean_sample[sampid].endswith('fq.gz'):
            code += 'zcat {cleandata} | sed -n \'1~4s/^/>/p;2~4p\' | \\\n'.format(cleandata = dict_clean_sample[sampid])
        elif dict_clean_sample[sampid].endswith('fa.gz'):
            code += 'zcat {cleandata} | \\\n'.format(cleandata = dict_clean_sample[sampid])
        else:
            code += 'cat {cleandata} | \\\n'.format(cleandata = dict_clean_sample[sampid])

        
        code +=( 
                '   awk \'NR%2==0{{if(length($0)>= {min_len} && length($0)<= {max_len}) print a"\\n"$0}}{{a=$0}}\' \\\n'
                '   >> {projpath}/1.QC/{sampid}/clean_data/{sampid}_remain_total.fa && \n'
                'sed -n \'2~2p\' {projpath}/1.QC/{sampid}/clean_data/{sampid}_remain_total.fa | \\\n'
                '   awk \'{{array[$0]++}}END{{for(key in array)print ">"key"("array[key]")""\\n"key}}\' \\\n'
                '   >> {projpath}/1.QC/{sampid}/clean_data/{sampid}_remain_uniq.fa && \n'
                '{perlExec} {BioScripts}/1QC/stat_seq_len_distribution.pl \\\n'
                '   -fa {sampid}_remain_uniq.fa \\\n'
                '   -pre {sampid}  -min {min_len}  -max {max_len} \\\n'
                '   -od  {projpath}/1.QC/{sampid}/clean_data && \n'
                'date +"%D %T -> Finsh Stat sRNA distribution"\n'
        ).format(perlExec = self.perlExec,BioScripts = self.BioScripts,min_len = self.min_len,
                 max_len = self.max_len,projpath = self.projpath,sampid = sampid)
       
        shell_path = os.path.join(self.projpath,'1.QC',sampid,'clean_data','{}_draw_length_distribution.sh'.format(sampid))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path, 1,'draw_length_distribution', sched=self.sched)
        self.job.add_atom_job(atom_job)

        
        
    def _category(self,sampid):
        code = (
            'date +"%D %T -> Start Category " && \\\n'
            '{perlExec} {BioScripts}/1QC/sRNA_compose.pl \\\n'
            '   -i {projpath}/1.QC/{sampid}/clean_data/{sampid}_remain_total.fa \\\n'
            '   -s {sampid} \\\n'
            '   -o {projpath}/1.QC/{sampid}/Category && \\\n'
            'date +"%D %T -> Finish Category "\n'
        ).format(projpath = self.projpath,sampid = sampid,perlExec = self.perlExec,
                 BioScripts = self.BioScripts)
        shell_path = os.path.join(self.projpath,'1.QC',sampid,'Category','{}_category.sh'.format(sampid))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path, 1,'category', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='{}_draw_length_distribution'.format(sampid))
        self.job.add_order(job_name,after_jobs='get_hairpin_mature_seq')


    def _MapStat(self,sampid):
        
        if 'tae' in self.abbr.lower():
            bowtieGenome = self.bowtie2+"/bowtie --large-index"                #小麦这种大基因组比对时候需要参加--large-index参数
        else:
            bowtieGenome = self.bowtie+"/bowtie"      

        code = (
            'date +"%D %T -> Stat Mapping and Stat " && \\\n\n'
            'if [ ! -f {projpath}/2.map/reference.mapping.stat ];then \n'
            '  echo -e  "Sample\\tTotal sRNA\\tMapped sRNA\\t\\"+\\"Mapped sRNA\\t\\"-\\"Mapped sRNA" \\\n'
            '  > {projpath}/2.map/reference.mapping.stat \n'
            'fi && \n\n'
            'awk \'{{if(/^>/){{split(substr($1,2),a,"(");split(a[2],b,")");'
            'id="{sampid}_"i++;for(j=1;j<=b[1];j++){{print ">"id"_x"b[1]"_"j"\\n"a[1]}}}}}}\' '
            '{projpath}/1.QC/{sampid}/clean_data/{sampid}_remain_uniq.fa > '
            '{projpath}/2.map/{sampid}/{sampid}.total.fa && \n\n'
            'awk \'{{if(/^>/){{split(substr($1,2),a,\"(\");split(a[2],b,\")\");'
            'id=\"{sampid}_\"i++;print \">\"id\"_x\"b[1]\"\\n\"a[1]}}}}\' '
            '{projpath}/1.QC/{sampid}/clean_data/{sampid}_remain_uniq.fa '
            '> {projpath}/2.map/{sampid}/{sampid}.collapse.fa && \n\n'
            '{bowtieGenome} -p 5 -v 1 -k 1 \\\n'
            '    {refer} \\\n'
            '    -f {projpath}/2.map/{sampid}/{sampid}.total.fa \\\n'
            '    --sam {projpath}/2.map/{sampid}/{sampid}.sam && \n\n'
            '{bowtieGenome} -p 5 -v 1 -k 1 \\\n'
            '    {refer} \\\n'
            '    -f {projpath}/2.map/{sampid}/{sampid}.total.fa \\\n'
            '    {projpath}/2.map/{sampid}/{sampid}.bwt \\\n'
            '    --un {projpath}/2.map/{sampid}/{sampid}.unmap.total.fa && \n\n'
            '{perlExec} {BioScripts}/2map/bwt2sta-uniq.pl  \\\n'
            '   {projpath}/2.map/{sampid}/{sampid}.bwt    \\\n'
            '   {projpath}/2.map/{sampid}/{sampid}.mapping.stat \\\n'
            '   {projpath}/2.map/{sampid}/{sampid}.map.collapse.fa && \n\n'
            'awk \'{{if(/^>/){{split($1,a,"_x");uniq++;count=a[2];total+=count}}else{{uniqbase+=length($1);'
            'totalbase+=count*length($1)}}}}END{{print "Total small RNA\\t"total"\\t"totalbase"\\t"uniq"\\t"uniqbase}}\' ' 
            '{projpath}/2.map/{sampid}/{sampid}.collapse.fa >>{projpath}/2.map/{sampid}/{sampid}.mapping.stat && \n\n'
            'awk -F"\\t" \'BEGIN{{sample="{sampid}"}}{{if(/^Total small RNA/){{total=$2}}else if(/^Total Mapped small RNA/)'
            '{{map=$2}}else if(/^Total Sense Mapped small RNA/){{sense=$2}}else if(/^Total Antisense Mapped small RNA/)'
            '{{antisense=$2}}}}END{{map_percent=sprintf("%.2f%",map/total*100);sense_percent=sprintf("%.2f%",sense/total*100);'
            'antisense_percent=sprintf("%.2f%",antisense/total*100);'
            'print sample"\\t"total" (100.00%)\\t"map" ("map_percent")\\t"sense" ("sense_percent")\\t"antisense" ("antisense_percent")"}}\' '
            '{projpath}/2.map/{sampid}/{sampid}.mapping.stat ' 
            '>>{projpath}/2.map/reference.mapping.stat && \n\n'
            'date +"%D %T -> Finish  Mapping and Stat "\n'
        ).format(projpath = self.projpath,sampid = sampid,perlExec = self.perlExec,
                 BioScripts = self.BioScripts,refer = self.refer,bowtie = self.bowtie,
                 bowtieGenome = bowtieGenome)

        shell_path = os.path.join(self.projpath,'2.map',sampid,'{}_MapStat.sh'.format(sampid))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path, 1,'MapStat', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='{}_draw_length_distribution'.format(sampid)) 
        self.job.add_order(job_name, after_jobs='get_hairpin_mature_seq')
        
        
    def _circos(self,sampid):
        code = (
                'date +"%D %T ->Start Genome distribution for {sampid}" && \\\n'
                'export PATH={circos}/bin:$PATH \n'
                'export PATH={R_v2153}:$PATH \n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:{LIBRARY_PATH}/HTSeq_Atlas:$LD_LIBRARY_PATH \n'
                '{samtools} view -bS {projpath}/2.map/{sampid}/{sampid}.sam  > {projpath}/2.map/{sampid}/{sampid}.bam && \n'
                '{python276} {BioScripts}/2map/circos_density_v2.3.py  \\\n'
                '   --bam {projpath}/2.map/{sampid}/{sampid}.bam \\\n'
                '   --fa {refer} \\\n'
                '   --n {chrNum} --iv 10000 --r 10000 \\\n'
                '   --outputfile {sampid}_circos \\\n'
                '   --outputdir {projpath}/2.map/{sampid}/Circos  && \n'
                'date +"%D %T ->Finish Genome distribution for {sampid}"\n'
        ).format(projpath = self.projpath,sampid = sampid,python276 = self.python276,
                 BioScripts = self.BioScripts,refer = self.refer,circos = self.circos,chrNum = self.chrNum,
                 R_v2153 = self.R_v2153, LIBRARY_PATH = self.LIBRARY_PATH,samtools = self.samtools)
        shell_path = os.path.join(self.projpath,'2.map',sampid,'Circos','{}_circos.sh'.format(sampid))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path, 1,'circos', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='{}_MapStat'.format(sampid))
#        self.job.add_order(job_name, after_jobs='srna_result') 
    
    def _MD5_Sum(self):
        pass
