#!/usr/bin/env python
#-*- coding:utf8 -*-
# Power by wangyunkai@novogene.com 2018-05-03 9:55:28

import os
import sys
import gzip
import glob
import shutil
import argparse

BASIC = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
        ))

sys.path.append(os.path.join(BASIC,'BioModule'))
from Stools import CheckTools
from Settings import BioConfig


class srnaResult(object):

    def __init__(self,**kwargs):
        self.projpath     = kwargs.get('projpath') if kwargs.get('projpath') else os.path.dirname(os.path.os.getcwd())
        self.project      = kwargs.get('project')
        self.org          = kwargs.get('org')
        self.samples      = kwargs.get('sample').strip().split(',')
        self.group        = CheckTools(kwargs).group
        self.groupname    = CheckTools(kwargs).groupname
        self.comparename  = CheckTools(kwargs).compare[1]
        self.venn_cluster = CheckTools(kwargs).venn_cluster[0]
        self.go           = kwargs.get('go')
        self.geneAnn      = kwargs.get('geneAnn')
        self.gene         = kwargs.get('gene')
        self.utr3         = kwargs.get('utr3')
        self.English      = kwargs.get('English')
        self.result       = os.path.join(self.projpath,self.project+"_sRNA_result",self.project+"_results")
        self._order       = 0
        self.BioScripts   = os.path.join(BASIC,"BioScripts")
        self.GOseq        = BioConfig.config.get('database','GOseq')
        self.perlExec     = BioConfig.config.get('srnaenv','perl_v5182')
        if self.English :
            self.readme   = os.path.join(self.BioScripts,"Readme","English")
        else:
            self.readme   = os.path.join(self.BioScripts,"Readme","Chinese")
        self.common       = kwargs.get('common')


    @property
    def order(self):
        self._order +=1
        return str(self._order)    
    

    def release(self):
        
        self.Example_data()
        self.qc()
        self.mapping()
        self.known_miRNA()
        self.ncRNA()
        #
        if self.org.strip().lower() == 'refplant':
            self.refplant()
        if 'norefplant' in self.org:
            self.norefplant()
        if self.org.strip().lower() == 'refanimal':
            self.refanimal()        
        if 'norefanimal' in self.org:
            self.norefanimal()
        #
        
        if len(self.samples) != 1:
            self.multi_sample_Enrichment()
        else:
            self.single_sample_Enrichment()
        #
        self.SuppFiles()


    def refplant(self):
        self.repeat()
        self.NAT()
        self.geneAnalysis()
        self.novel_miRNA()
        self.TAS()
        self.Category()
        self.miRNA_editing()
        self.miRNA_family()
        self.diff_analy()
        self.miRNA_target_plant()


    def norefplant(self):
        self.novel_miRNA()
        self.Category()
        self.miRNA_editing()
        self.miRNA_family()
        self.diff_analy()
        self.miRNA_target_plant()

    def refanimal(self):
        self.repeat()
        self.geneAnalysis()
        self.novel_miRNA()
        self.Category()
        self.miRNA_editing()
        self.miRNA_family()
        self.diff_analy()
        self.miRNA_target_refanimal()

    def norefanimal(self): 
        self.novel_miRNA()
        self.Category()
        self.miRNA_editing()
        self.miRNA_family()
        self.diff_analy()
        self.miRNA_target_norefanimal()        


    def SuppFiles(self):

        SuppFilesdir = self.make_dir('SuppFiles')
        self.copy_dir(os.path.join(
                      self.projpath,'11.edit_family','edit_analy','known_novel_hairpin.fa'
                                  ),os.path.join(SuppFilesdir,'hairpin.fa')
                     )
        assert not os.system('cat {known_mature} {novel_mature} > {total_mature}'.format(
                                known_mature = os.path.join(self.result,'*Known_miRNA','mature.fa'),
                                novel_mature = os.path.join(self.result,'*Novel_miRNA','mature.fa'),
                                total_mature = os.path.join(SuppFilesdir,'mature.fa'))
                             )
                                     
        assert not os.system('{perlExec} {BioScripts}/Release/GO_class.pl {go} {GOseq} {outfile}'.format(
                                 perlExec = self.perlExec, BioScripts = self.BioScripts,go = self.go,
                                 GOseq = self.GOseq,outfile = os.path.join(SuppFilesdir,'gene.goannot.xls'))
                           )
                
        self.copy_dir(os.path.join(
                      self.projpath,'*.diff','diffAnalysisResult','{meanscount.txt,meanstpm.txt}'
                                  ),SuppFilesdir
                     )
        self.copy_dir(self.geneAnn,os.path.join(SuppFilesdir,'gene.description.xls'))
        self.copy_dir(os.path.join(
                      self.readme,'0.SuppFiles.Readme.txt'
                                  ),SuppFilesdir
                     )
        
        assert not os.system(
                             '{perlExec} {BioScripts}/Release/known_haipin_mature_result.pl \\\n'
                             '  -readcount {readcount} \\\n'
                             '  -h {hairpin} \\\n'
                             '  -m {mature} \\\n'
                             '  -pairs {hairpin_mature_pairs} \\\n'
                             '  -o {outfile} '.format(
                                      perlExec  = self.perlExec,BioScripts = self.BioScripts,
                                      readcount = os.path.join(self.projpath,'3.known/known_miRNAs/mature.readcount'),
                                      hairpin   = os.path.join(self.projpath,'3.known/known_miRNAs/hairpin.fa'),
                                      mature    = os.path.join(self.projpath,'3.known/known_miRNAs/mature.fa'),
                                      hairpin_mature_pairs = os.path.join(self.projpath,'3.known/known_miRNAs/hairpin_mature.pairs'),
                                      outfile   = os.path.join(SuppFilesdir,'known_hairpin_mature.xls')
                                                      )
                            )

        assert not os.system(
                             '{perlExec} {BioScripts}/Release/novel_haipin_mature_result.pl \\\n'
                             '  -readcount {readcount} \\\n'
                             '  -h {hairpin} \\\n'
                             '  -m {mature} \\\n'
                             '  -str {hairpin_str} \\\n'
                             '  -pos {hairpin_pos} \\\n'
                             '  -o {outfile} '.format(
                                      perlExec  = self.perlExec,BioScripts = self.BioScripts,
                                      readcount = os.path.join(self.projpath,'8.novel/novel_miRNAs/mature.readcount'),
                                      hairpin   = os.path.join(self.projpath,'8.novel/novel_miRNAs/hairpin.fa'),
                                      mature    = os.path.join(self.projpath,'8.novel/novel_miRNAs/mature.fa'),
                                      hairpin_str = os.path.join(self.projpath,'8.novel/novel_miRNAs/hairpin.str'),
                                      hairpin_pos = os.path.join(self.projpath,'8.novel/novel_miRNAs/hairpin.pos'),
                                      outfile   = os.path.join(SuppFilesdir,'novel_hairpin_mature.xls')
                                                     )
                            )

        if 'refanimal' in self.org:
            assert not os.system('{perlExec} {BioScripts}/Release/targets_pairs.transfer.v3.pl {inf} {ouf}'.format(
                                         perlExec  = self.perlExec,BioScripts = self.BioScripts,
                                         inf = os.path.join(self.result,'*miRNA_target/all_target_gene.xls'),
                                         ouf = os.path.join(SuppFilesdir,'targets.pairs.allcom')
                                         )
                                )
        else:
            assert not os.system('{perlExec} {BioScripts}/Release/targets_pairs.transfer.v3.pl {inf} {ouf}'.format(
                                         perlExec  = self.perlExec,BioScripts = self.BioScripts,
                                         inf = os.path.join(self.result,'*miRNA_target/*_targets.pairs'),
                                         ouf = os.path.join(SuppFilesdir,'targets.pairs.allcom')
                                         )
                                )
            
        assert not os.system("sed -e 's/,/\\t/g' {inf} > {ouf}".format(
                                      inf = os.path.join(SuppFilesdir,'targets.pairs.allcom'),
                                      ouf = os.path.join(SuppFilesdir,'targets.pairs.alltab')
                                      )
                            )
        
        if 'noref' not in self.org:
            if 'refplant' in self.org:
                self.copy_dir(self.gene,os.path.join(SuppFilesdir,'gene.fa'))
            else:
                self.copy_dir(self.utr3,os.path.join(SuppFilesdir,'3_utr.fa'))

                
    def Example_data(self):
        
        #RawDatadir,CleanDatadir,order = self.make_dir('Example_data',['1.RawData','2.CleanData'])
        CleanDatadir,order = self.make_dir('Example_data',['1.CleanData'])       

        ''' 
        for eachsample in self.samples:
            files = glob.glob(self.projpath+'/raw_data/'+eachsample+'/'+eachsample+'.fq.gz')[0]
            f_in = gzip.open(files,'rb')
            f_out = [f_in.readline() for num in range(20)]
            f_out_name = os.path.join(RawDatadir,eachsample+'.example.raw.fastq')
            open(f_out_name,'w').writelines(f_out)
            f_in.close()
        shutil.copy('{}/1.1RawData.txt'.format(self.readme),os.path.join(RawDatadir,order+'.1.RawData.Readme.txt'))
        '''        

        for eachsample in self.samples:
            files=glob.glob(self.projpath+'/1.QC/'+eachsample+'/clean_data/'+eachsample+'_remain_total.fa')[0]
            f_in=open(files,'rb')
            f_out=[f_in.readline() for num in range(20)]
            f_out_name = os.path.join(CleanDatadir,eachsample+'.example.clean.fastq')
            open(f_out_name,'w').writelines(f_out)
            f_in.close()
        #shutil.copy('{}/1.2CleanData.txt'.format(self.readme),os.path.join(CleanDatadir,order+'.2.CleanData.Readme.txt'))
        shutil.copy('{}/1.2CleanData.txt'.format(self.readme),os.path.join(CleanDatadir,order+'.1.CleanData.Readme.txt'))


    def qc(self):

        #RawData_ErrorRate,RawData_Stat,ReadsClassification,Length_Filter,order = self.make_dir(
        #    'QualityControl',['1.RawData_ErrorRate','2.RawData_Stat','3.ReadsClassification','4.Length_Filter']
        #)
        
        Length_Filter,order = self.make_dir('QualityControl',['1.Length_Filter'])

        #qc_report = os.path.join(self.projpath,self.project+'_QC_results','results')
        #dirname = os.path.dirname(RawData_ErrorRate)
        #self.copy_dir(os.path.join(qc_report,'1RawData_ErrorRate/*'),RawData_ErrorRate)
        #self.copy_dir(os.path.join(qc_report,'2RawData_Stat/*'),RawData_Stat)
        #self.copy_dir(os.path.join(qc_report,'3ReadsClassification/*'),ReadsClassification)
        self.copy_dir(os.path.join(self.projpath,'1.QC/*/clean_data/*len_distribution.p*'),Length_Filter)
        #shutil.copy(
        #            '{}/2.1RawData_ErrorRate.txt'.format(self.readme),
        #            os.path.join(RawData_ErrorRate,order+'.1.RawData_ErrorRate.txt')
        #            )
        #shutil.copy(
        #            '{}/2.2RawData_Stat.txt'.format(self.readme),
        #            os.path.join(RawData_Stat,order+'.2.RawData_Stat.txt')
        #            )
        #shutil.copy(
        #            '{}/2.3ReadsClassification.txt'.format(self.readme),
        #            os.path.join(ReadsClassification,order+'.3.ReadsClassification.txt')
        #            )
        shutil.copy(
                    '{}/2.4Length_Filter.txt'.format(self.readme),
                    os.path.join(Length_Filter,order+'.1.Length_Filter.txt')
                    )
        #shutil.copy(
        #            '{}/2QC_2.txt'.format(self.readme),
        #            os.path.join(dirname,order+'.QC.Readme.txt')
        #            )     
        if 'y' in self.common :
            assert not os.system('mkdir -p {}'.format(os.path.join(dirname,order+'.5.Common_Specific_sRNA')))
            self.copy_dir(os.path.join(self.projpath,'1.QC','*_vs_*/*venn*'),dirname)
            shutil.copy(
                        '{}/2.5Common-Specific_sRNA.txt'.format(self.readme),
                        os.path.join(dirname,order+'.5.Common_Specific_sRNA',order+'.5.Common_Specific_sRNA.Readme.txt')
                        )

    def mapping(self):

        Mapping_Stat,order = self.make_dir('Mapping_Stat')
        analydir = os.path.join(self.projpath,'2.map')
        shutil.copy(
                    os.path.join(analydir,'reference.mapping.stat'),
                    Mapping_Stat
                   )
        for eachsample in self.samples:
            shutil.copy(
                        os.path.join(analydir,eachsample,eachsample+'.mapping.stat'),
                        Mapping_Stat
                       )
            if 'noref' not in self.org:
                shutil.copy(
                            os.path.join(analydir,eachsample,'Circos',eachsample+'_circos.png'),
                            os.path.join(Mapping_Stat,eachsample+'.circos.png')
                           )
                shutil.copy(
                            os.path.join(analydir,eachsample,'Circos',eachsample+'_circos.svg'),
                            os.path.join(Mapping_Stat,eachsample+'.circos.svg')
                           )        
        self.copy_dir(os.path.join(
                                   self.readme,'3Mapping_Stat.txt'
                                  ),os.path.join(Mapping_Stat,order+'.Mapping_Stat.Readme.txt')
                     )


    def known_miRNA(self):
        
        known_miRNA,order = self.make_dir('Known_miRNA')
        analydir = os.path.join(self.projpath,'3.known')        
        assert not os.system('mkdir -p {}'.format(
                             os.path.join(known_miRNA,'Structure_plot_example'))
                            )
        assert not os.system('tar -zcf {tardir} -C {rawdir} {file}'.format(
                             tardir = os.path.join(known_miRNA,'Structure_plot.tar.gz'),
                             rawdir = os.path.join(analydir,'known_miRNAs'),
                             file = 'image')
                             )
        assert not os.system('ls {} 2>&1 | head -n 10 | xargs -n1 -i cp {{}} {}'.format(
                            os.path.join(analydir,'known_miRNAs','image/*jpg'),
                            os.path.join(known_miRNA,'Structure_plot_example'))
                            )

        self.copy_dir(os.path.join(
                      analydir,'known_miRNAs','{hairpin*,mature*,known*.map.*,miRBase.mrd}'
                                  ),known_miRNA
                     )
                             
        for eachsample in self.samples:        
            shutil.copy(
                        os.path.join(analydir,eachsample+'.firstbase.png'),
                        known_miRNA
                       )
            shutil.copy(
                        os.path.join(analydir,eachsample+'.position.png'),
                        known_miRNA
                       )
        self.copy_dir(os.path.join(
                                   self.readme,'4Known_miRNA.txt'
                                  ),os.path.join(known_miRNA,order+'.Known_miRNA.Readme.txt')
                     )

        
    def ncRNA(self):

        ncRNA,order = self.make_dir('ncRNA')
        analydir = os.path.join(self.projpath,'4.ncRNA')
        self.copy_dir(os.path.join(
                      analydir,'output','{*.map.fas,rc.stat,uc.stat}'
                                  ),ncRNA
                     )
        self.copy_dir(os.path.join(
                      self.readme,'5ncRNA.txt'
                                  ),os.path.join(ncRNA,order+'.ncRNA.Readme.txt')
                     )                    

    def repeat(self):
        repeat,order = self.make_dir('Repeat')
        analydir = os.path.join(self.projpath,'5.repeat')
        self.copy_dir(os.path.join(
                      analydir,'{repeat.map.fas,repeat.rc.stat,repeat.uc.stat,*bar.png,*bar.pdf}'
                                  ),repeat
                     )
        self.copy_dir(os.path.join(
                      self.readme,'6Repeats.txt'
                                  ),os.path.join(repeat,order+'.Repeats.Readme.txt')
                     )

    def NAT(self):
        NAT,order = self.make_dir('NAT')
        analydir = os.path.join(self.projpath,'6.NAT')
        self.copy_dir(os.path.join(
                      analydir,'output','{*.map.fas,rc.stat,uc.stat,NAT.rc.stat,NAT.uc.stat}'
                                  ),NAT
                     )
        self.copy_dir(os.path.join(
                      self.readme,'7NAT.txt'
                                  ),os.path.join(NAT,order+'.NAT.Readme.txt')
                     )

    def geneAnalysis(self):
        gene,order = self.make_dir('gene')
        analydir = os.path.join(self.projpath,'7.gene')
        self.copy_dir(os.path.join(
                      analydir,'output','{*.map.fas,rc.stat,uc.stat}'
                                  ),gene
                     )
        self.copy_dir(os.path.join(
                      self.readme,'8gene.txt'
                                  ),os.path.join(gene,order+'.gene.Readme.txt')
                     )

    def novel_miRNA(self):

        Novel_miRNA,order = self.make_dir('Novel_miRNA')
        analydir = os.path.join(self.projpath,'8.novel')        
        assert not os.system('mkdir -p {}'.format(
                             os.path.join(Novel_miRNA,'Structure_plot_example'))
                            )
        assert not os.system('tar -zcf {tardir} -C {rawdir} {file}'.format(
                             tardir = os.path.join(Novel_miRNA,'Structure_plot.tar.gz'),
                             rawdir = os.path.join(analydir,'novel_miRNAs'),
                             file = 'image')
                             )
        assert not os.system('ls {} 2>&1 | head -n 10 | xargs -n1 -i cp {{}} {}'.format(
                            os.path.join(analydir,'novel_miRNAs','image/*jpg'),
                            os.path.join(Novel_miRNA,'Structure_plot_example'))
                            )
        self.copy_dir(os.path.join(
                      analydir,'novel_miRNAs','{hairpin*,mature*,novel*.map.*,miRBase.mrd,star*}'
                                  ),Novel_miRNA
                     )
                             
        for eachsample in self.samples:        
            shutil.copy(
                        os.path.join(analydir,eachsample+'.firstbase.png'),
                        Novel_miRNA
                       )
            shutil.copy(
                        os.path.join(analydir,eachsample+'.position.png'),
                        Novel_miRNA
                       )
        self.copy_dir(os.path.join(
                                   self.readme,'9Novel_miRNA.txt'
                                  ),os.path.join(Novel_miRNA,order+'.Novel_miRNA.Readme.txt')
                     )
            
    def TAS(self):
        
        tas,order = self.make_dir('TAS')
        analydir = os.path.join(self.projpath,'9.TAS')
        self.copy_dir(os.path.join(
                      analydir,'output','{TAS.map.fas,TAS.rc.stat,TAS.uc.stat}'
                                  ),tas
                     )
        self.copy_dir(os.path.join(
                      self.readme,'10TAS.txt'
                                  ),os.path.join(tas,order+'.TAS.Readme.txt')
                     )
        
    def Category(self):

        category,order = self.make_dir('Category')
        analydir = os.path.join(self.projpath,'10.Category')
        self.copy_dir(os.path.join(
                      analydir,'{*_full.txt,*pie.png,*pie.pdf}'
                                  ),category
                     )
        self.copy_dir(os.path.join(
                      self.readme,'11Category.txt'
                                  ),os.path.join(category,order+'.Category.Readme.txt')
                     )

    def miRNA_editing(self):

        miRNA_editing,order = self.make_dir('miRNA_editing')
        analydir = os.path.join(self.projpath,'11.edit_family')
        self.copy_dir(os.path.join(
                      analydir,'edit_analy','*editing_stats*.txt'
                                  ),miRNA_editing
                     )
        self.copy_dir(os.path.join(
                      self.readme,'12miRNA_editing.txt'
                                  ),os.path.join(miRNA_editing,order+'.miRNA_editing.Readme.txt')
                     )
        
    
    def miRNA_family(self):

        miRNA_family,order = self.make_dir('miRNA_family')
        analydir = os.path.join(self.projpath,'11.edit_family')
        self.copy_dir(os.path.join(
                      analydir,'family_analy','*_family.*.txt'
                                  ),miRNA_family
                     )
        self.copy_dir(os.path.join(
                      self.readme,'13miRNA_family.txt'
                                  ),os.path.join(miRNA_family,order+'.miRNA_family.Readme.txt')
                     )


    def diff_analy(self):

        analydir = os.path.join(self.projpath,'13.diff')
        if len(self.samples) == 1:
            miRNAExp,miRNAExpdensity,order = self.make_dir('DiffExprAnalysis',
                    ['1.miRNAExp','2.miRNAExpdensity'])
        else:
            if self.venn_cluster:
                (miRNAExp,miRNAExpdensity,CorAnalysis,DiffExprAnalysis,
                     DEsFilter,DEcluster,DEvenn,order) = self.make_dir('DiffExprAnalysis',[
                        '1.miRNAExp','2.miRNAExpdensity','3.CorAnalysis','4.DiffExprAnalysis',
                        '5.DEsFilter','6.DEcluster','7.DEvenn'
                         ])
            else:
                (miRNAExp,miRNAExpdensity,CorAnalysis,DiffExprAnalysis,
                     DEsFilter,DEcluster,order) = self.make_dir('DiffExprAnalysis',[
                        '1.miRNAExp','2.miRNAExpdensity','3.CorAnalysis',
                        '4.DiffExprAnalysis','5.DEsFilter','6.DEcluster'
                         ])

        self.copy_dir(os.path.join(
                        analydir,'diffAnalysisResult','{Readcount_TPM.xls,TPM_interval.xls}'
                            ),miRNAExp
                     )
        self.copy_dir(os.path.join(
                      self.readme,'14.1miRNAExp.txt'
                                  ),os.path.join(miRNAExp,order+'.1.miRNAExp.Readme.txt')
                     )
        
        self.copy_dir(os.path.join(
                        analydir,'diffAnalysisResult','{TPM_density_distribution.*,TPM_boxplot.*}'
                            ),miRNAExpdensity
                     )
        self.copy_dir(os.path.join(
                      self.readme,'14.2miRNAExpdensity.txt'
                                  ),os.path.join(miRNAExpdensity,order+'.2.miRNAExpdensity.Readme.txt')
                     )
        
        if len(self.samples) != 1:
        
            self.copy_dir(os.path.join(
                            analydir,'diffAnalysisResult','corr_plot/*'
                                ),CorAnalysis
                         )
            self.copy_dir(os.path.join(
                          self.readme,'14.3CorAnalysis.txt'
                                      ),os.path.join(CorAnalysis,order+'.3.CorAnalysis.Readme.txt')
                         )
            
            for each_compare in self.comparename.replace(':','vs').split(','):
                assert not os.system('mkdir -p {}'.format(os.path.join(DiffExprAnalysis,each_compare))) 
                self.copy_dir(os.path.join(
                                analydir,'diffAnalysisResult',each_compare,'{*vs*xls,*vs*txt}'
                                          ),os.path.join(DiffExprAnalysis,each_compare)
                             )
                self.copy_dir(os.path.join(
                                analydir,'diffAnalysisResult',each_compare,'*Volcanoplot*'
                                          ),DEsFilter
                             )                

            self.copy_dir(os.path.join(
                          self.readme,'14.4DiffExprAnalysis.txt'
                                      ),os.path.join(DiffExprAnalysis,order+'.4.DiffExprAnalysis.Readme.txt')
                         )
        
            self.copy_dir(os.path.join(
                          self.readme,'14.5DEsFilter.txt'
                                      ),os.path.join(DEsFilter,order+'.5.DEsFilter.Readme.txt')
                         )
            
            #self.copy_dir(os.path.join(
            #                analydir,'diffAnalysisResult',
            #                '{*heatmap.png,*heatmap.pdf,DE_union_for_cluster,SOM_cluster,K_means_cluster}'
            #                    ),DEcluster
            #             )
            #self.copy_dir(os.path.join(
            #              self.readme,'14.6DEcluster.txt'
            #                          ),os.path.join(DEcluster,order+'.6.DEcluster.Readme.txt')
            #             )

            if self.venn_cluster:
                
                self.copy_dir(os.path.join(
                                analydir,'diffAnalysisResult','venn'
                                    ),DEvenn
                             )
                assert not os.system('rm -rf {}'.format(os.path.join(
                                    DEvenn,'venn','{type.txt,Signature.xls,plot_venn.R}'))
                                )
                self.copy_dir(os.path.join(
                              self.readme,'14.7DEvenn.txt'
                                          ),os.path.join(DEvenn,order+'.7.DEvenn.Readme.txt')
                             )

            if len(self.groupname.split(','))==2:
                
                tmpdir = os.path.dirname(os.path.abspath(miRNAExp))
                DEvenn = os.path.join(tmpdir,order+'.7.DEvenn')
                assert not os.system('mkdir -p {}'.format(DEvenn))
                               
                self.copy_dir(os.path.join(
                                analydir,'diffAnalysisResult','venn'
                                    ),DEvenn
                             )
                self.copy_dir(os.path.join(
                              self.readme,'14.8DEvenn.txt'
                                          ),os.path.join(DEvenn,order+'.7.DEvenn.Readme.txt')                             )
                
                
    def miRNA_target_refanimal(self):
        
        miRNA_target,order = self.make_dir('miRNA_target')
        analydir = os.path.join(self.projpath,'12.target')
        self.copy_dir(os.path.join(
                      analydir,'Common','{commom_target.xls,all_target_gene.xls,commom_target_example.xls,all_target_gene.xls.annotate}'
                                  ),miRNA_target
                     )

        assert not os.system('gzip -c {}/miRanda/miranda_targets_out.fmt > {}/miranda_targets_out.fmt.gz'.format(analydir,miRNA_target))
        assert not os.system('gzip -c {}/RNAhybrid/RNAhybrid_miRNA_target_pairs > {}/RNAhybrid_miRNA_target_pairs.gz'.format(analydir,miRNA_target))       
        assert not os.system('gzip -c {}/PITA/PITA_pita_results_targets.tab > {}/PITA_pita_results_targets.tab.gz'.format(analydir,miRNA_target)) 

        self.copy_dir(os.path.join(
                      self.readme,'15miRNA_target.txt'
                                  ),os.path.join(miRNA_target,order+'.miRNA_target.Readme.txt')
                     )

    def miRNA_target_plant(self):
        
        miRNA_target,order = self.make_dir('miRNA_target')
        analydir = os.path.join(self.projpath,'12.target')
        self.copy_dir(os.path.join(
                    analydir,'{*targets.txt,*_example.xls,*targets.pairs.annotate,*targets.pairs}'
                            ),miRNA_target
                     )
        self.copy_dir(os.path.join(
                      self.readme,'15miRNA_target.txt'
                                  ),os.path.join(miRNA_target,order+'.miRNA_target.Readme.txt')
                     )  
    
        
    def miRNA_target_norefanimal(self):
        
        miRNA_target,order = self.make_dir('miRNA_target')
        analydir = os.path.join(self.projpath,'12.target')
        self.copy_dir(os.path.join(
                    analydir,'{*targets,*_example.xls,*targets.pairs.annotate,*targets.pairs}'
                            ),miRNA_target
                     )
        assert not os.system('gzip -c {}/mature.miranda_targets_out.fmt > {}/mature.miranda_targets_out.fmt.gz'.format(analydir,miRNA_target))
        self.copy_dir(os.path.join(
                    self.readme,'15miRNA_target.txt'
                            ),os.path.join(miRNA_target,order+'.miRNA_target.Readme.txt')
                     )


    def single_sample_Enrichment(self):
        
        Enrichment,order = self.make_dir('Enrichment')
        analydir = os.path.join(self.projpath,'14.enrich')
        godir    = os.path.join(Enrichment,self.samples[0],'GOenrichment') 
        keggdir  = os.path.join(Enrichment,self.samples[0],'KEGGenrichment')
        assert not os.system('mkdir -p {}'.format(godir))
        assert not os.system('mkdir -p {}'.format(keggdir))
        self.copy_dir(os.path.join(
                      analydir,'GO_enrichment','{*DEG*,*png,*pdf}'
                                  ),godir
                     )
        self.copy_dir(os.path.join(
                      analydir,'Pathway','{*scatterplot.png,*scatterplot.pdf,src}'
                                  ),keggdir
                     )
        self.copy_dir(os.path.join(
                      analydir,self.samples[0],'Pathway','*html'),
                      os.path.join(
                      keggdir,self.samples[0]+'.DEG_enriched_KEGG_pathway_API.html')
                      )
        self.copy_dir(os.path.join(
                      analydir,'Pathway',self.samples[0]+'.identify.xls'),
                      os.path.join(
                      keggdir,self.samples[0]+'.DEG_KEGG_pathway_enrichment_result.xls')
                      )
        self.copy_dir(os.path.join(
                      analydir,'Pathway','top*.xls'),
                      os.path.join(
                      keggdir,self.samples[0]+'.DEG_enriched_KEGG_pathway_top20.xls')
                     )
        self.copy_dir(os.path.join(
                      self.readme,'16Enrichment.txt'
                                  ),os.path.join(os.path.dirname(GOenrichment),order+'.Enrichment.Readme.txt')
                     )

    def multi_sample_Enrichment(self):

        Enrichment,order = self.make_dir('Enrichment')
        analydir = os.path.join(self.projpath,'14.enrich')
        self.copy_dir(os.path.join(
                      self.readme,'16Enrichment.txt'
                                  ),os.path.join(Enrichment,order+'.Enrichment.Readme.txt')
                     )

        for each_compare in self.comparename.replace(':','vs').split(','):
            temdir = os.path.join(Enrichment,each_compare)
            assert not os.system('mkdir -p {}'.format(temdir))
            self.copy_dir(os.path.join(
                          analydir,each_compare,'{*diffmiRNA-geneid,*diffmiRNAID,*diffmiRNA-gene.pairs}'
                                      ),temdir
                         )
            tmpgodir   = os.path.join(temdir,'GOenrichment')
            assert not os.system('mkdir -p {}'.format(tmpgodir))
            self.copy_dir(os.path.join(
                          analydir,each_compare,'GO2','*DEG*'
                                      ),tmpgodir
                          )
            os.system('rm -rf {}'.format(os.path.join(tmpgodir,
                                  '{*classification2_gene_count.txt,*classification2.pdf,*classification2.png}')
                                        )
                     )

            tmpkeggdir = os.path.join(temdir,'KEGGenrichment')
            assert not os.system('mkdir -p {}'.format(tmpkeggdir))
            self.copy_dir(os.path.join(
                          analydir,each_compare,'Pathway','{*scatterplot.pdf,*scatterplot.png,src}'
                                      ),tmpkeggdir
                         )
            
            self.copy_dir(os.path.join(
                          analydir,each_compare,'Pathway','*html'),
                          os.path.join(
                          tmpkeggdir,each_compare+'.DEG_enriched_KEGG_pathway_API.html')
                         )

            if 'noref' not in self.org:
 
                self.copy_dir(os.path.join(
                              analydir,each_compare,'Pathway',each_compare+'.identify.xls'),
                              os.path.join(
                              tmpkeggdir,each_compare+'.DEG_KEGG_pathway_enrichment_result.xls')
                              )
                self.copy_dir(os.path.join(
                              analydir,each_compare,'Pathway','top*.xls'),
                              os.path.join(
                              tmpkeggdir,each_compare+'.DEG_enriched_KEGG_pathway_top20.xls')
                              )
                
            else:
                
                self.copy_dir(os.path.join(
                              analydir,each_compare,'Pathway',each_compare+'*result.xls'),
                              os.path.join(
                              tmpkeggdir,each_compare+'.DEG_KEGG_pathway_enrichment_result.xls')
                              )
                self.copy_dir(os.path.join(
                              analydir,each_compare,'Pathway','*_top20.xls'),
                              os.path.join(
                              tmpkeggdir,each_compare+'.DEG_enriched_KEGG_pathway_top20.xls')
                              )


    def make_dir(self,analysis,subdirs=None):

        subdirL = list()
        order = self.order
        analydir = os.path.join(self.result,order + '.' + analysis)
        if 'SuppFiles' in analysis:
            SupFdir = os.path.join(self.result,'0.SuppFiles')
            assert not os.system(
                                'mkdir -p {}'.format(SupFdir)
                                )
            return SupFdir
        else:
            if subdirs:
                for subdir in subdirs:
                    absubdir = os.path.join(analydir,order + '.' + subdir)
                    assert not os.system('mkdir -p {}'.format(absubdir))
                    subdirL.append(absubdir)
                subdirL.append(order)
                return subdirL
            else:
                assert not os.system('mkdir -p {}'.format(analydir))
                return analydir,order

    @staticmethod
    def copy_dir(rawdir,newdir):
        if os.path.isfile(rawdir):
            assert not os.system(
                        'cp -f {rawdir} {newdir}'.format(
                        rawdir = rawdir, newdir = newdir)
                        )
        else:
            ## assert not os.system(
            os.system(
                        'cp -rf {rawdir} {newdir}'.format(
                        rawdir = rawdir, newdir = newdir)
                        )


def Argument():
    parser = argparse.ArgumentParser(
    description = 'srna pipline result v1.0',
    prog='srnaResult',
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
    parser.add_argument(
        '--go', 
        help='go file, used for "14.enrich" analysis', 
        default=True
    )
    parser.add_argument(
        '--geneAnn', 
        help='gene anno file, need title line, used for "12.target" analysis', 
        default=None
    )
    parser.add_argument(
        '--gene',
        help='gene.fa, used for "ref 12.plant target" analysis', 
        default=None
    )
    parser.add_argument(
        '--utr3', 
        help='3UTR.fa, used for "ref 12.animal target" analysis', 
        default=None
    )
    argv = vars(parser.parse_args())
    return argv


def main():
    args = Argument() 
    rad = srnaResult(**args)    
    rad.release()

if __name__ == '__main__':
    main()
