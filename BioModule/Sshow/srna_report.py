#!/usr/bin/env python
#-*- coding:utf8 -*-
# Power by wangyunkai@novogene.com 2018-05-03 14:53:21

import re
import os
import sys
import time
import glob
import shutil
import linecache
import argparse
from collections import defaultdict

reload(sys)
sys.setdefaultencoding('utf-8')
from django.template import Template, Context, loader
from django.conf import settings


BASIC = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
        ))
sys.path.append(os.path.join(BASIC,'BioModule'))
from Stools import CheckTools
from Settings import BioConfig


class srnaReport(object):
        
    def __init__(self,**kwargs):
        self.projpath      = kwargs.get('projpath') if kwargs.get('projpath') else os.path.dirname(os.path.os.getcwd())
        self.project       = kwargs.get('project')
        self.org           = kwargs.get('org')
        self.samples       = kwargs.get('sample').strip().split(',')
        self.BioScripts    = BioConfig.config.get('srnaenv','BioScripts')
        self.common        = kwargs.get('common')
        self.group         = CheckTools(kwargs).group
        self.groupname     = CheckTools(kwargs).groupname
        self.compare       = CheckTools(kwargs).compare[0]
        self.comparename   = CheckTools(kwargs).compare[1]
        self.venn_cluster  = CheckTools(kwargs).venn_cluster[0]
        self.venn_cluster_name = CheckTools(kwargs).venn_cluster[1]
        self.result        = os.path.join(self.projpath,self.project+"_sRNA_result",self.project+"_results")
        self.report        = os.path.join(self.projpath,self.project+"_sRNA_result",self.project+"_report")
        self.contract      = kwargs.get('contract')
        if not os.path.exists(self.report):
            os.makedirs(self.report)
        self.render        = defaultdict(dict)
        self.src           = os.path.join(self.report,'src')
        self.perlExec      = BioConfig.config.get('srnaenv','perl_v5182')
        self.English       = kwargs.get('English')
        self._code_number  = kwargs.get('project')
        self._html_order   = 0
        self.wkhtmltopdf   = BioConfig.config.get('software','wkhtmltopdf')


    _TUE  = '\033[36m'
    _FLE  = '\033[31m'
    _ENDC = '\033[0m'
    _WARNING = '\033[33m'


    @property
    def code_number(self):
        pattern = re.compile(r'.*-[BS].-\d+')
        if pattern.findall(self.project):
            self._code_number = pattern.findall(self.project)[0]
        else:
            #print self._FLE,"please check your project Parameter ,follow the rule : 项目编号-分析号-部门编号-物种名称",self._ENDC
            self._code_number = self.project
        return self._code_number

    @property
    def html_order(self):
        self._html_order +=1
        return str(self._html_order)


    def run(self):
        print  self.project + ' report generating...'
        self.reportTemplate()
        self.reportRender()
        self.qc()
        self.mapping()
        self.known_miRNA()
        self.analyncRNA()
        if 'noref' in self.org:
            self.noref()
        else:
            if 'refplant' in self.org:
                self.refplant()
            elif 'refanimal' in self.org:
                self.refanimal()
            else:
                self.noref()

        self.Category()
        self.analyedit()
        self.analyfamily()
        self.analydif()
        self.analyTarget()
        self.analyEnrich()
        self.getReport()
        #self.htmltopdf()

    def refplant(self):
        self.analyrepeat()
        self.analyNAT()
        self.analyGene()
        self.novel_miRNA()
        self.analyTAS()

    def refanimal(self):
        self.analyrepeat()
        self.analyGene()
        self.novel_miRNA()

    def noref(self):
        self.novel_miRNA()
        
    def reportTemplate(self):
        os.chdir(self.result)
        if self.English:
            settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, TEMPLATE_DIRS=(
                                    '{}/BioReport/English/Report_template_2.0/'.format(self.BioScripts),)
                               )
            assert not os.system('cp -r {BioScripts}/BioReport/English/Report_template_2.0/src {report}'.format(
                                            BioScripts = self.BioScripts,report = self.report)
                                )
            assert not  os.system('tree -d -v -L 2 -o {report}/DirectoryTree.html \\'
                                  '-H ../../../{project}_results'.format(report = self.report,
                                                                         project = self.project)
                                 )
            assert not  os.system('{perlExec} {BioScripts}/BioReport/result_tree_EN.pl \\\n'
                                  '   {report}/DirectoryTree.html \\\n'
                                  '   {report}/src/html/DirectoryTree.html'.format(
                                      perlExec = self.perlExec,BioScripts = self.BioScripts,report = self.report)
                                 )
            assert not  os.system('rm -rf {}/DirectoryTree.html'.format(self.report))
        else: 
            settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, TEMPLATE_DIRS=(
                                    '{}/BioReport/Chinese/Report_template_2.2/'.format(self.BioScripts),)
                              )
            assert not os.system('cp -r {BioScripts}/BioReport/Chinese/Report_template_2.2/src {report}'.format(
                                            BioScripts = self.BioScripts,report = self.report)
                                )
            #assert not os.system('cp -r {BioScripts}/BioReport/Novofinder_sRNA_manual.pdf {src}/pdf'.format(
            #                                    BioScripts = self.BioScripts, src = self.src)
            #                    )
            #assert not  os.system('tree -d -v -L 2 -o {report}/DirectoryTree.html \\'
            #                      '-H ../../../{project}_results'.format(report = self.report,
            #                                                             project = self.project)
            #                     )
            #assert not  os.system('{perlExec} {BioScripts}/BioReport/result_tree.pl \\\n'
            #                      '   {report}/DirectoryTree.html \\\n'
            #                      '   {report}/src/html/DirectoryTree.html'.format(
            #                          perlExec = self.perlExec,BioScripts = self.BioScripts,report = self.report)
            #                     )
            #assert not  os.system('rm -rf {}/DirectoryTree.html'.format(self.report))
        os.chdir(self.report)


    def getReport(self):

    
        t=loader.get_template('yankegen_smallrna_report.templete.html')
        c=Context(self.render)
        html=t.render(c)
        open(os.path.join(self.report,'yankegen_smallrna_report.html'),'w').write(html)
 
        '''
        t=loader.get_template('index.html')
        c=Context(self.render)
        html=t.render(c)
        open(os.path.join(self.report,self.code_number+'_sRNA_report.html'),'w').write(html)

        t=loader.get_template('src/html/right.html')
        c=Context(self.render)
        html=t.render(c)
        open('{src}/html/right.html'.format(src = self.src),'w').write(html)

        t=loader.get_template('src/html/top.html')
        c=Context(self.render)
        html=t.render(c)
        open('{src}/html/top.html'.format(src = self.src),'w').write(html) 
        '''
        


    def reportRender(self):
                
        self.render['title']         = self.project
        self.render['time']          = time.strftime('%Y-%m-%d',time.localtime(time.time()))
        contract_list                = self.contract.strip().split('_')
        contract_id                  = contract_list[0]
        contract_name                = contract_list[1]
        self.render['contract_id']   = contract_id
        self.render['contract_name'] = contract_name        
        self.render['code_number']   = self.code_number
        #summary

        known_miRNA_stat = glob.glob(os.path.join(self.result,'*Known_miRNA','known_miRNA.map.stat'))[0] 
        self.render['table_known_map_summary_head'] = self.head2html(linecache.getlines(known_miRNA_stat)[0])
        self.render['table_known_map_summary'] = self.tab2html(linecache.getlines(known_miRNA_stat)[1:3])
        novel_miRNA_stat = glob.glob(os.path.join(self.result,'*Novel_miRNA','novel_miRNA.map.stat'))[0]
        self.render['table_novel_map_summary_head'] = self.head2html(linecache.getlines(novel_miRNA_stat)[0])
        self.render['table_novel_map_summary'] = self.tab2html(linecache.getlines(novel_miRNA_stat)[1:3])
        if len(self.samples) !=1 and self.compare:
            self.render['flag_diff'] = True 
            table_diffsum = glob.glob(os.path.join(self.projpath,'*.diff','diff_sum.txt'))[0]
            self.render['table_diffsum_head'] = self.head2html(linecache.getlines(table_diffsum)[0])
            self.render['table_diffsum']      = self.tab2html(linecache.getlines(table_diffsum)[1:])


    def qc(self):

        raw_order = self.html_order
        self.render['html_raw_order']   = raw_order
        qc_order  = self.html_order       
        self.render['html_qc_order']    = qc_order
        self.render['figure_qc_order']  = qc_order
        #
        '''
        self.renderFigure('figure_error','2.QualityControl/2.1.RawData_ErrorRate',
                            '_error_rate_distribution.png','.error_rate_distribution.png')
        '''
        self.renderFigure('figure_lenFilter','2.QualityControl/2.1.Length_Filter',
                             '_seq_len_distribution.png','_seq_len_distribution.png')
        #
        '''
        self.render['table_qc_order']   = qc_order
        table_basic                     = os.path.join(self.result,'2.QualityControl','2.2.RawData_Stat','RawData_Stat.xls')
        self.render['table_basic']      = self.tab2array(linecache.getlines(table_basic)[1:])
        #
        table_filter                    = os.path.join(self.result,'2.QualityControl','2.3.ReadsClassification','clean_process_overview.xls')
        self.render['table_filter']     = self.tab2array(linecache.getlines(table_filter)[1:])
        #
        '''
        table_lenFilter                 = os.path.join(self.result,'2.QualityControl','2.1.Length_Filter','total_uniq.xls')
        self.render['table_lenFilter']  = self.tab2array(linecache.getlines(table_lenFilter)[1:])
        #
        if 'y' in self.common:
            self.render['flag_common']  = True
            self.render['figure_seqVenn_total'] = list()
            self.render['figure_seqVenn_uniq']  = list()
            for group in [ sampid+"_vs_"+sample_list[index+1] 
                           for index,sampid in list(enumerate(sample_list))[:-1]
                         ]: 
                k=5
                assert not os.system('convert -resize 600 {Common_Specific_sRNA}/{group}_venn_chart_Tol.png '
                                     '{report}/src/images/{group}_total_venn.png'.format(
                                           RawData_ErrorRate = os.path.join(self.result,'2.QualityControl','2.5.Common_Specific_sRNA'),
                                           eachsample = eachsample, report = self.report)
                                    )
                self.render['figure_seqVenn_total'].append(["'"+'src/images/'+each+'_total_venn.png'+"'"])

                assert not os.system('convert -resize 600 {Common_Specific_sRNA}/{group}_venn_chart_uni.png '
                                     '{report}/src/images/{group}_uniq_venn.png'.format(
                                           RawData_ErrorRate = os.path.join(self.result,'2.QualityControl','2.5.Common_Specific_sRNA'),
                                           eachsample = eachsample, report = self.report)
                                    )
                self.render['figure_seqVenn_uniq'].append(["'"+'src/images/'+each+'_uniq_venn.png'+"'"])
                k-=1
                if k==0:
                    break


    def mapping(self):
        
        self.render['flag_map']=True
        map_order = self.html_order
        self.render['html_map_order']   = map_order
        self.render['table_map_order']  = map_order
        #
        table_map = glob.glob(os.path.join(self.result,'*Mapping_Stat','reference.mapping.stat'))[0]
        self.render['table_map']   = self.tab2array(linecache.getlines(table_map)[1:])
        #
        if 'noref' in self.org:
            self.render['known_circos'] = False
        else:
            if 'refanimal' in self.org:
                self.render['flag_animal']  = True
            elif 'refplant' in self.org:
                self.render['flag_plant']   = True
            self.render['flag_ref']         = True
            self.render['figure_map_order'] = map_order
            self.render['known_circos']     = True
            #
            self.renderFigure('figure_map','*Mapping_Stat','.circos.png','.circos.png')


    def known_miRNA(self):
        
        self.render['flag_known'] = True
        known_order = self.html_order
        self.render['html_known_order'] = known_order
        self.render['table_known_miRNA_order'] = known_order
        #
        table_known_map = glob.glob(os.path.join(self.result,'*Known_miRNA','known_miRNA.map.stat'))[0]
        self.render['table_known_map_head'] = self.head2html(linecache.getlines(table_known_map)[0])
        self.render['table_known_map']      = self.tab2html(linecache.getlines(table_known_map))
        #
        self.render['figure_known_order']   = known_order
        self.render['figure_known']         = ['src/images/known_miRNA_structure.jpg']
        #
        table_known_rc = glob.glob(os.path.join(self.result,'*Known_miRNA','mature.readcount'))[0]
        self.render['table_known_rc_head']  = self.head2html(linecache.getlines(table_known_rc)[0])
        self.render['table_known_rc']       = self.tab2html(linecache.getlines(table_known_rc)[1:],5)
        #
        known_mrd = glob.glob(os.path.join(self.result,'*Known_miRNA','miRBase.mrd'))[0]
        self.render['known_mrd']            = self.tab2html(linecache.getlines(known_mrd),10)
        #
        self.renderFigure('figure_known_bias1','*Known_miRNA','.firstbase.png','.known_firstbase.png')
        self.renderFigure('figure_known_bias2','*Known_miRNA','.position.png','.known_position.png')
   
 
    def analyncRNA(self):
        
        self.render['flag_ncRNA'] = True
        ncRNA_order = self.html_order
        self.render['html_ncRNA_order']   = ncRNA_order
        self.render['table_ncRNA_order']  = ncRNA_order
        #
        table_ncRNA = glob.glob(os.path.join(self.result,'*ncRNA','rc.stat'))[0]
        self.render['table_ncRNA_head']  = self.head2html(linecache.getlines(table_ncRNA)[0])
        self.render['table_ncRNA']       = self.tab2html(linecache.getlines(table_ncRNA)[1:])


    def analyrepeat(self):
        
        self.render['flag_repAna'] = True
        repAna_order = self.html_order
        self.render['html_repAna_order']   = repAna_order 
        self.render['figure_repAna_order'] = repAna_order
        #
        self.renderFigure('figure_repAna1','*.Repeat','_bar.png','_total_repeats.png','rc.stat_')
        self.renderFigure('figure_repAna2','*.Repeat','_bar.png','_uniq_repeats.png','uc.stat_')


    def analyNAT(self):

        self.render['flag_NAT'] = True
        NAT_order = self.html_order
        self.render['html_NAT_order']   = NAT_order
        self.render['table_NAT_order']  = NAT_order
        #
        table_NAT = glob.glob(os.path.join(self.result,'*.NAT','rc.stat'))[0]
        self.render['table_NAT_head']  = self.head2html(linecache.getlines(table_NAT)[0])
        self.render['table_NAT']       = self.tab2html(linecache.getlines(table_NAT)[1:])
        

    def analyGene(self):
        
        self.render['flag_gene'] = True
        gene_order = self.html_order
        self.render['html_gene_order']  = gene_order
        self.render['table_gene_order'] = gene_order
        #
        table_gene = glob.glob(os.path.join(self.result,'*.gene','rc.stat'))[0]
        self.render['table_gene_head']  = self.head2html(linecache.getlines(table_gene)[0])
        self.render['table_gene']       = self.tab2html(linecache.getlines(table_gene)[1:])


    def novel_miRNA(self):
        
        self.render['flag_novel'] = True
        novel_order = self.html_order
        self.render['html_novel_order'] = novel_order
        self.render['table_novel_miRNA_order'] = novel_order
        #
        table_novel_map = glob.glob(os.path.join(self.result,'*.Novel_miRNA','novel_miRNA.map.stat'))[0]
        self.render['table_novel_map_head']  = self.head2html(linecache.getlines(table_novel_map)[0])
        self.render['table_novel_map']       = self.tab2html(linecache.getlines(table_novel_map)[1:])
        #
        self.render['figure_novel_miRNA_order'] = novel_order
        self.render['figure_novel'] = ['src/images/novel_miRNA_structure.jpg']
        #
        table_novel_rc = glob.glob(os.path.join(self.result,'*.Novel_miRNA','mature.readcount'))[0]
        self.render['table_novel_rc_head']  = self.head2html(linecache.getlines(table_novel_rc)[0])
        self.render['table_novel_rc']       = self.tab2html(linecache.getlines(table_novel_rc)[1:],5)
        #
        novel_mrd = glob.glob(os.path.join(self.result,'*.Novel_miRNA','miRBase.mrd'))[0]
        self.render['novel_mrd']             = self.tab2html(linecache.getlines(novel_mrd),10)
        #
        self.renderFigure('figure_novel_bias1','*.Novel_miRNA','.firstbase.png','_novel_firstbase.png')
        self.renderFigure('figure_novel_bias2','*.Novel_miRNA','.position.png','_novel_position.png')


    def analyTAS(self):
        
        self.render['flag_TAS'] = True
        TAS_order  = self.html_order
        self.render['html_TAS_order']  = TAS_order
        self.render['table_TAS_order'] = TAS_order
        #
        table_TAS = glob.glob(os.path.join(self.result,'*.TAS','TAS.rc.stat'))[0]
        self.render['table_TAS_head']  = self.head2html(linecache.getlines(table_TAS)[0])
        self.render['table_TAS']       = self.tab2html(linecache.getlines(table_TAS)[1:])

    def Category(self):
        
        self.render['flag_category'] = True
        category_order = self.html_order
        self.render['html_category_order']  = category_order
        self.render['table_category_order'] = category_order
        #
        table_category = glob.glob(os.path.join(self.result,'*.Category','category_rc_full.txt'))[0]
        self.render['table_category_head']  = self.head2html(linecache.getlines(table_category)[0])
        self.render['table_category']       = self.tab2html(linecache.getlines(table_category)[1:])


    def analyedit(self):
        
        self.render['flag_edit'] = True
        edit_order = self.html_order
        self.render['html_edit_order'] = edit_order
        #
        table_edit = glob.glob(os.path.join(self.result,'*.miRNA_editing','*editing_stats.example.txt'))[0]
        self.render['edit'] = self.tab2html(linecache.getlines(table_edit))


    def analyfamily(self):
        
        self.render['flag_family']=True
        family_order = self.html_order
        self.render['html_family_order']  = family_order
        self.render['table_family_order'] = family_order
        #
        table_family = glob.glob(os.path.join(self.result,'*.miRNA_family','*_miRNA_family.mir_sign.example.txt'))[0]
        self.render['table_family_head']  = 'Species'+self.head2html(linecache.getlines(table_family)[0])
        self.render['table_family'] = self.tab2html(linecache.getlines(table_family)[1:])

    def analydif(self):
        
        diff_order = self.html_order
        self.render['html_diff_order']   = diff_order
        self.render['table_diff_order']  = diff_order
        #
        table_tpm = glob.glob(os.path.join(self.result,'*.DiffExprAnalysis','*.miRNAExp','Readcount_TPM.xls'))[0]
        self.render['table_tpm_head']  = self.head2html(linecache.getlines(table_tpm)[0])
        self.render['table_tpm']       = self.tab2html(linecache.getlines(table_tpm)[1:],10)
        #
        self.render['figure_diff_order'] = diff_order
        self.render['figure_tpm']  = list()
        self.renderSingleFigure('figure_tpm','*.DiffExprAnalysis/*.miRNAExpdensity','TPM_density_distribution.png','TPM_density_distribution.png')
        #
        if self.render['flag_diff']:
            self.render['figure_cor']  = list()
            self.renderSingleFigure('figure_cor','*.DiffExprAnalysis/*.CorAnalysis','cor_pearson.png','cor_pearson.png')
            k=5
            for png in glob.iglob(self.result+'/*.DiffExprAnalysis/*.CorAnalysis/*.scatter.png'):
                name = re.search(r'(/.*/)(.*scatter)\.png',png).group(2)
                assert not os.system('convert -resize 600 %s %s' % (png,self.report+'/src/images/'+name+'.png'))
                self.render['figure_cor'].append(["'"+'src/images/'+os.path.basename(png)+"'"])
                k-=1
                if k==0:
                    break
            #
            compare_exp = self.comparename.split(',')[0].replace(':','vs')
            table_diff = glob.glob(os.path.join(self.result,'*.DiffExprAnalysis/*.DiffExprAnalysis',compare_exp,compare_exp+'.DE.xls'))[0]
            self.render['table_diff_head']  = self.head2html(linecache.getlines(table_diff)[0])
            self.render['table_diff']       = self.tab2html(linecache.getlines(table_diff)[1:],5)      
            #
            self.render['figure_volcano'] = list()
            for eachcom in self.comparename.split(','):
                temp = eachcom.replace(':','vs')
                self.renderSingleFigure('figure_volcano','*.DiffExprAnalysis/*.DEsFilter',temp+'.Volcanoplot.png',temp+'.Volcanoplot.png')
            #
            self.render['figure_heatmap'] = list()
            self.renderSingleFigure('figure_heatmap','*.DiffExprAnalysis/*.DEcluster','Hcluster_heatmap.png','Hcluster.png')
            #
            venn_cluster_vs_names = list()
            if self.venn_cluster_name != None:
                for each in self.venn_cluster_name.split(','):
                    if each.count(':') >= 2 and each.count(':') <= 5:
                        self.render['flag_venn'] = True
                        venn_cluster_vs_names.append(each.replace(':','vs'))
            #
            self.render['figure_venn'] = list()
            if self.render['flag_venn']:
                for each in venn_cluster_vs_names:
                    self.renderSingleFigure('figure_venn','*.DiffExprAnalysis/*.DEvenn/venn/'+each,each+'.venn.png',each+'.DEG_Venn_diagram.png')
            #
            if len(self.groupname.split(',')) == 2:
                groupnames = self.groupname.split(',')
                self.render['flag_venn'] = True
                self.renderSingleFigure('figure_venn','*.DiffExprAnalysis/*.DEvenn/venn/'+groupnames[0]+'_'+groupnames[1],
                                               groupnames[0]+'_'+groupnames[1]+'.DEG_Venn_diagram.png',
                                               groupnames[0]+'_'+groupnames[1]+'.DEG_Venn_diagram.png')

    def analyTarget(self):
         
        self.render['flag_target'] = True
        target_order = self.html_order
        self.render['html_target_order'] = target_order
        #
        table_target = glob.glob(os.path.join(self.result,'*.miRNA_target','*_example.xls'))[0]
        self.render['table_target_head']  = self.head2html(linecache.getlines(table_target)[0])
        self.render['table_target']       = self.tab2html(linecache.getlines(table_target)[1:],10)

    def analyEnrich(self):
        
        self.render['flag_enrich'] = True
        enrich_order = self.html_order
        self.render['html_enrich_order']  = enrich_order
        self.render['table_enrich_order'] = enrich_order
        #
        if len(self.samples) == 1:
            self.render['flag_single'] = True
            table_go = glob.glob(os.path.join(self.result,'*.Enrichment/*/GOenrichment','*.DEG_GO_enrichment_result.xls'))[0] 
        else:
            table_go = glob.glob(os.path.join(self.result,'*.Enrichment/*vs*/GOenrichment','*vs*.DEG_GO_enrichment_result.xls'))[0]
        self.render['table_go_head']  = self.head2html('\t'.join(linecache.getlines(table_go)[0].split('\t')))
        self.render['table_go']       = self.tab2html(self.extractColumn(linecache.getlines(table_go)[1:]),5)
        #
        if self.render['flag_single']:
            self.render['figure_enrich_order'] = enrich_order
            #
            self.render['figure_goBar'] = list()            
            self.renderSingleFigure('figure_goBar','*.Enrichment/*/GOenrichment',
                                        '*.DEG_Enriched_GO_classification.png',self.samples[0]+'.Enriched_GO_classification.png')
            #
            self.render['figure_goDAG'] = list()
            self.renderSingleFigure('figure_goDAG','*.Enrichment/*/GOenrichment',
                                        '*.DEG_Enriched_GO_cc_DAG.png',self.samples[0]+'.Enriched_GO_cc_DAG.png')
            self.renderSingleFigure('figure_goDAG','*.Enrichment/*/GOenrichment'
                                        '*.DEG_Enriched_GO_bp_DAG.png',self.samples[0]+'.Enriched_GO_bp_DAG.png')
            self.renderSingleFigure('figure_goDAG','*.Enrichment/*/GOenrichment'
                                        '*.DEG_Enriched_GO_mf_DAG.png',self.samples[0]+'.Enriched_GO_mf_DAG.png')
            #
            enrich_file = os.path.join(self.result,'*.Enrichment/*/KEGGenrichment','*.DEG_KEGG_pathway_enrichment_result.xls')
            flag = 0
            k = 5
            sig_map = ''
            kegg_png = list()
            self.render['table_kegg'] = list()
            for eachLine in open(glob.glob(enrich_file)[0]):
                if eachLine.startswith('#Term'):
                    flag = 1
                    continue
                if flag == 1 and eachLine.strip() != '' and not eachLine.startswith('--'):
                    temper = eachLine.split('\t')
                    if sig_map == '':
                        sig_map = temper[2].strip()
                        sig_map_png = glob.glob(os.path.join(self.result,'*.Enrichment/*/KEGGenrichment/src/',sig_map+'.png'))[0]
                        if os.path.isfile(sig_map_png):
                            kegg_png.append(sig_map_png)
                            sig_map = ''
                        else:
                            sig_map = ''
                    self.render['table_kegg'].append(temper[:7])
                    k-=1
                if k == 0:
                    break            
            #
            self.render['figure_kegg_scat'] = list()
            self.renderSingleFigure('figure_kegg_scat','*.Enrichment/*/KEGGenrichment',
                                      '*.DEG_enriched_KEGG_pathway_scatterplot.png',self.samples[0]+'.KEGG_pathway_scatterplot.png')
            self.render['figure_kegg_path'] = list()
            for png in kegg_png:
                name=re.search(r'(/.*/)(.*)\.png',png).group(2)
                assert not os.system('convert -resize 600 %s %s' % (png,self.report+'/src/images/'+name+'.png'))
                self.render['figure_kegg_path'].append(["'"+'src/images/'+os.path.basename(png)+"'"])
        else:
            self.render['figure_enrich_order'] = enrich_order
            self.render['figure_goBar'] = list()
            self.render['figure_goDAG'] = list() 
            for eachcom in self.comparename.split(','):
                temp = eachcom.replace(':','vs')
                self.renderSingleFigure('figure_goBar','*.Enrichment/'+temp+'/GOenrichment',
                                         temp+'.DEG_Enriched_GO_classification.png',temp+'.Enriched_GO_classification.png')
                #
                self.renderSingleFigure('figure_goDAG','*.Enrichment/'+temp+'/GOenrichment',
                                            temp+'.DEG_Enriched_GO_cc_DAG.png',temp+'.Enriched_GO_cc_DAG.png')
                self.renderSingleFigure('figure_goDAG','*.Enrichment/'+temp+'/GOenrichment',
                                            temp+'.DEG_Enriched_GO_bp_DAG.png',temp+'.Enriched_GO_bp_DAG.png')
                self.renderSingleFigure('figure_goDAG','*.Enrichment/'+temp+'/GOenrichment',
                                            temp+'.DEG_Enriched_GO_mf_DAG.png',temp+'.Enriched_GO_mf_DAG.png')
            #
            flag=0
            sig_map = ''
            kegg_png = list()
            self.render['table_kegg'] = list()
            for eachcom in self.comparename.split(','):
                temp = eachcom.replace(':','vs')
                kegg_png_tmp = glob.glob(self.result+'/*.Enrichment/'+temp+'/KEGGenrichment/src/*.png')
                if kegg_png_tmp:
                    enrich_compare= temp
                    break
                else:
                    continue
            if len(kegg_png_tmp)>5:
                k=5
            else:
                k=len(kegg_png_tmp)
            #
            if enrich_compare:
                enrich_file = self.result+'/*.Enrichment/'+enrich_compare+'/KEGGenrichment/*result.xls'
            for eachLine in open(glob.glob(enrich_file)[0]):
                if eachLine.startswith('#Term'):
                    flag=1
                    continue
                elif flag == 1 and eachLine.strip() != '' and not eachLine.startswith('-'):
                    temper=eachLine.split('\t')
                    if sig_map == '':
                        sig_map = temper[2].strip()
                        sig_map_png = glob.glob(self.result+'/*.Enrichment/*vs*/KEGGenrichment/src/'+sig_map+'.png')[0]
                        if os.path.isfile(sig_map_png):
                            kegg_png.append(sig_map_png)
                            sig_map = ''
                        else:
                            sig_map = ''
                    self.render['table_kegg'].append(temper[:7])
                    k-=1
                if k == 0:
                    break
            #
            self.render['figure_kegg_scat'] = list()
            for eachcom in self.comparename.split(','):
                temp = eachcom.replace(':','vs')
                self.renderSingleFigure('figure_kegg_scat','*.Enrichment/'+temp+'/KEGGenrichment',
                                            '*.DEG_enriched_KEGG_pathway_scatterplot.png',temp+'.KEGG_pathway_scatterplot.png')
            #
            self.render['figure_kegg_path'] = list()
            for png in kegg_png:
                name=re.search(r'(/.*/)(.*)\.png',png).group(2)
                assert not os.system('convert -resize 600 %s %s' % (png,self.report+'/src/images/'+name+'.png'))
                self.render['figure_kegg_path'].append(["'"+'src/images/'+os.path.basename(png)+"'"])


    def htmltopdf(self):
        
        code = ('{wkhtmltopdf} --page-width 350mm --page-height 495mm \\'
                '-n --print-media-type --footer-center \'[page] / [topage]\' \\' 
                'cover {src}/html/cover.html toc --toc-header-text \'content\' \\'
                '--toc-text-size-shrink 1 {src}/html/right.html {report}/{code_number}_sRNA_report.pdf').format(
                     wkhtmltopdf = self.wkhtmltopdf,src = self.src,
                     report = self.report,code_number = self.code_number)

        assert not os.system(code)             
        assert not os.system("sed -i 's/#00//g' {}".format(os.path.join(self.report,self.code_number+'_sRNA_report.pdf')))


    @staticmethod
    def extractColumn(lines):
        nline = list()
        for line in lines:
            nline.append('\t'.join(line.split('\t')[:9]))
        return nline

    @staticmethod
    def head2html(head):
        headHtml='<tr><th>'+'</th><th>'.join(head.strip().split('\t'))+'</th></tr>'
        return headHtml

    @staticmethod
    def tab2html(lines,num=None):
        tab2list = list()
        if num and len(lines)>=num:
            for line in lines[:num]:
                tabHtml='<tr><td>'+'</td><td>'.join(line.strip().split('\t'))+'</td></tr>'
                tab2list.append(tabHtml)
        else:
            for line in lines:
                tabHtml='<tr><td>'+'</td><td>'.join(line.strip().split('\t'))+'</td></tr>'
                tab2list.append(tabHtml)

        return tab2list
   
    @staticmethod
    def tab2array(lines):
        tab2array = list()
        for line in lines:
            tab2array.append(line.strip().split('\t'))
        return tab2array                

    def renderFigure(self,figureName,analysis,PicAnalySuffix,PicReportSuffix,PicAnalyPrefix=None):
        
        self.render[figureName] = list()
        PicAnalyPrefix = PicAnalyPrefix if PicAnalyPrefix else '' 
        for eachsample in self.samples:
            os.system('convert -resize 600 {analydir}/{PicPreName} '
                        '{report}/src/images/{PicNowName}'.format(
                             analydir   = os.path.join(self.result,analysis),
                             PicPreName = PicAnalyPrefix + eachsample + PicAnalySuffix,
                             PicNowName = eachsample + PicReportSuffix,
                             report = self.report)
                     )
            self.render[figureName].append(["'"+'src/images/'+eachsample+PicReportSuffix+"'"])

    def renderSingleFigure(self,figureName,analysis,PicAnalySuffix,PicReportSuffix,PicAnalyPrefix=None):
        PicAnalyPrefix = PicAnalyPrefix if PicAnalyPrefix else ''
        os.system('convert -resize 600 {analydir}/{PicPreName} '
                    '{report}/src/images/{PicNowName}'.format(
                         analydir   = os.path.join(self.result,analysis),
                         PicPreName = PicAnalyPrefix + PicAnalySuffix,
                         PicNowName = PicReportSuffix,
                         report = self.report)
                 )
        self.render[figureName].append(["'"+'src/images/'+PicReportSuffix+"'"])


def Argument():
    parser = argparse.ArgumentParser(
    description = 'srna pipline report v1.0',
    prog='srnaReport',
    formatter_class=argparse.RawTextHelpFormatter,
    epilog='    Contact:  rag@novogene.com'
    ) 
    parser.add_argument(
        '--projpath', 
        metavar = 'path', help="The project path", 
        default=None
    )
    parser.add_argument(
        '--project',
        help=('project name, maybe same with the name of root dir,\n'
              'which will be displayed in the final report title'),
        required=True
    )
    parser.add_argument(
        '--sample',
        help="sample names(sample1,sample2,...)",
        required=True
    )
    parser.add_argument(
        '--English',
        help='English report or readme',
        default=None
    )
    parser.add_argument(
        '--org',
        help='organism',
        required=True
    )
    parser.add_argument(
        '--common',
        help='common specific',
        default='n',
        required=None
    )
    parser.add_argument(
        '--contract',
        help='contract name',
        required=True
    )
    parser.add_argument(
        '--group',
        help=('sample grouping way, e.g. TR1,TR2,TS1,TS2,TR1:TR2,TS1:TS2,\n'
              '\':\' used in intra-group;  \',\'used in inter-group'),
        default=None
    )
    parser.add_argument(
        '--groupname',
        help='group name, splited by "," , e.g. TR1,TR2,TS1,TS2,TR,TS',
        default=None
    )
    parser.add_argument(
        '--compare',
        help=('group compare way, e.g. \'2:1,1:3,2:3\'\n'
              '1,2,3 were groupname order,1:2(treat:control),\n' 
              'intra-group splited by \':\',inter-group by \',\''),
        default=None
    )
    parser.add_argument(
        '--venn_cluster',
        help=('venn ploting way; suit for 2~4 compare groups; \n'
              '\'_\' split compare group in the same one plot; \n'
              '\',\' split different plots, e.g. \'2:1_1:3_2:3,1:3_2:3\'\n'),
        default=None
    )    
    argv = vars(parser.parse_args())
    return argv

def main():
    args = Argument()   
    rad = srnaReport(**args)
    rad.run()
    
if __name__ == '__main__':
    main()

