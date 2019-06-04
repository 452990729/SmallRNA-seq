#!/usr/bin/env python
# _*_ coding:utf-8 _*_

import os
import sys
import argparse
from datetime import datetime

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR+'/BioModule')

from Settings import BioConfig
from QcMapping  import QcMapping
from Stools import Job,AtomJob
from sRNAClassify import sRNAClassify
from miRNAanaly import miRNAanaly


class run_pipeline(object):
    
    def __init__(self,args):
        self.args      =  args
        self.projpath  =  args.get('projpath')
        self.new_job   =  args.get('new_job')
        if not os.path.exists(os.path.join(args.get('projpath'),'log')):
            os.makedirs(os.path.join(args.get('projpath'),'log'))
        self.log       =  os.path.join(args.get('projpath'),'log')
        self.job       =  Job(os.path.join(self.log,self.new_job))

    def Run(self):
        job_file = os.path.join(self.projpath,self.new_job)
        QcMapping(self.args,self.job).run()
        sRNAClassify(self.args,self.job).run()
        miRNAanaly(self.args,self.job).run()
        self.job.write(job_file) 


def Argument():
    
    parser = argparse.ArgumentParser(
        description = 'smallrna pipline Reformation from RAD',
        prog='smallrna',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Contact:  wangyunkai@novogene.com'
    )
    parser.add_argument(
        '--projpath',
        metavar='Str',
        help='The path of analysis directory[defalut="."]',
        default=os.getcwd()
    )
    parser.add_argument(
        '--project',
        help=('project name, maybe same with the name of root dir \n' 
              'which will be displayed in the final report title'),
        required=True
    )
    parser.add_argument(
        '--org',
        help='organism',
        required=True
    )
    parser.add_argument(
        '--contract',
        help='contract name',
        required=True
    )
    parser.add_argument(
        '--mapfile',
        help='mapfile files 2column: "CleanDataFataGZfile\\tSampleName"',
        required=True
    )
    parser.add_argument(
        '--out_dir',
        help="the outputdir of analysis result",
        default=None
    )
    parser.add_argument(
        '--ad5',
        default='GTTCAGAGTTCTACAGTCCGACGATC',
        help="5' adapter (default=%(default)s)",
    ) 
    parser.add_argument(
        '--ad3',
        default='AGATCGGAAGAGCACACGTCT',
        help="3' adapter (default=%(default)s)", 
    )
    parser.add_argument(
        '--sample',
        help="sample names(sample1,sample2,...)",
        required=True
    )
    parser.add_argument(
        '--group',
        help=("sample grouping way,e.g.TR1,TR2,TS1,TS2,TR1:TR2,TS1:TS2 \n"
              "':' used in intra-group;','used in inter-group"),
        default=None
    )
    parser.add_argument(
        '--groupname',
        help='group name, splited by "," ,e.g. TR1,TR2,TS1,TS2,TR,TS',
        default=None
    )
    parser.add_argument(
        '--compare',
        help=("group compare way, e.g.'2:1,1:3,2:3' 1,2,3 were groupname order \n"
              "1:2 (treat:control),intra-group splited by ':',inter-group by ','"),
        default=None
    )
    parser.add_argument(
        '--venn_cluster',
        help=("venn ploting way; suit for 2~4 compare groups;\n" 
              "'_' split compare group in the same one plot;\n"
              "',' split different plots, e.g. '2:1_1:3_2:3,1:3_2:3'"),
        default=None
    )
    parser.add_argument(
        '--mdspe',
        help='mode species name',
        default=None
    )
    parser.add_argument(
        '--refer', 
        help='reference for mapping, may be genome or assembled transcript', 
        default=None
    )
    parser.add_argument(
        '--chrNum', 
        help=('chr number to plot, recommend 10 when chr is contig or scafford \n'
              'used for "2.map" analysis'),
        default=10
    )
    parser.add_argument(
        '--abbr', 
        help='abbreviation of species from miRBase', 
        required=True
    )
    parser.add_argument(
        '--rRNA', 
        help='rRNA.fa', 
        default=None
    )
    parser.add_argument(
        '--tRNA', 
        help='tRNA.fa', 
        default=None
    )
    parser.add_argument(
        '--snRNA', 
        help='snRNA.fa', 
        default=None
    )
    parser.add_argument(
        '--snoRNA', 
        help='snoRNA.fa', 
        default=None
    )
    parser.add_argument(
        '--spe', 
        help=('species or related species full name, selected from \n' 
              '(BioModule/5repeat/RepeatMasker_Species.txt)'), 
        default=None
    )
    parser.add_argument(
        '--repdir', 
        help='repeat prediction dir', 
        default=None
    ) 
    parser.add_argument(
        '--NAT', 
        help=('NAT species abbrevation, selected from (PlantNATsDB/spe.list),'
              'if not exsit, use "other" as substitute'), 
        default=None
    )
    parser.add_argument(
        '--exon', 
        help='exon.fa, used for "7.gene" analysis', 
        default=None
    )
    parser.add_argument(
        '--intron', 
        help='intron.fa, used for "7.gene" analysis', 
        default=None
    )
    parser.add_argument(
        '--mode', 
        help='mode for predict novel miRNA, 1 for animal; 2 for monocot; 3 for dicots', 
        choices=['1','2','3'],
        default='1'
    ) 
    parser.add_argument(
        '--gene',
        help='gene.fa, used for "ref 9.plant target" analysis', 
        default=None
    )
    parser.add_argument(
        '--utr3', 
        help='3UTR.fa, used for "ref 9.animal target" analysis', 
        default=None
    )
    parser.add_argument(
        '--gff3', 
        help='CDS gff3 file, used for "noRef 9.animal target" analysis', 
        default=None
    )
    parser.add_argument(
        '--geneAnn', 
        help='gene anno file, need title line, used for "9.target" analysis', 
        default=None
    )
    parser.add_argument(
        '--gtf', 
        help='gtf file, used for "10.enrich" analysis', 
        default=None
    )
    parser.add_argument(
        '--go', 
        help='go file, used for "10.enrich" analysis', 
        default=True
    )
    parser.add_argument(
        '--kegg', 
        help='KEGG species abbrevation, used for "10.enrich" analysis', 
        default=None
    )
    parser.add_argument(
        '--ko', 
        help='ko file, used for "10.enrich" analysis', 
        default=None
    )
    parser.add_argument(
        '--length', 
        help=('length file, used for "10.enrich" analysis,'
              'which is GeneINFO in trinity dir'), 
        default=None
    )
    parser.add_argument(
        '--English',
        help='English report',
        default=None
    )
    parser.add_argument(
        '--ownername',
        help='xinxi name',
        required=None
    )
    parser.add_argument(
        '--yunying',
        help='yunying de name',
        required=None
    )
    parser.add_argument(
        '--common',
        help='common specific',
        default='n',
        required=None
    )
    parser.add_argument(
        '--exosome',
        help='if exosome project,y,else n',
        default='n'
    )
    parser.add_argument(
        '--type',
        help=('type of species,include 3utr_human,3utr_fly,3utr_worm,'
              'for animal to predict target gene'),
        default=None
    )
    parser.add_argument(
        '--analy_code',
        metavar='Str',
        help='The anaylsis code list seperated by comma',
        default=None
    )
    parser.add_argument(
        '--new_job', 
        default="%s.job" % datetime.strftime(datetime.now(), '%Y.%m.%d.%H'),
        help='use for construct job name and log dir, default is year.month.date.hour.job'
    )
    parser.add_argument(
        "--sched", 
        default=None,
        help='sjm投递的参数，默认为根据qselect提取用户可用队列构造,为了防止跟命令行参数冲突，本参数传递时必须加引号\n'
        '例如："-V -cwd -S /bin/bash -q rad.q "'
    )
    parser.add_argument(
        "--mdspedb",
        default='{mdspedb}'.format(mdspedb = BioConfig.config.get('database','ModeSpeDB')),
        help='model species database storage path (default=%(default)s)'
    )
    argv = vars(parser.parse_args())
    return argv


def main():
    args = Argument()
    run_pipeline(args).Run()


if __name__ == '__main__':
    main()
