#!/usr/bin/env python
# _*_ coding:utf-8 _*_
 

import os
import glob
from Settings import BioConfig
from Stools import CheckTools, AtomJob, Job, write_shell ,make_dir
from collections import defaultdict


class miRNAanaly(object):

    def __init__(self,args,job):

        self.args         = args
        self.job          = job
        self.projpath     = args.get('projpath') or os.getcwd()
        self.mirdeep2     = BioConfig.config.get('software','mirdeep2')
        self.sample_list  = CheckTools(args).samples
        self.sample       = args.get('sample')
#        self.group        = CheckTools(args).group
#        self.groupname    = CheckTools(args).groupname
#        self.compare      = CheckTools(args).compare[0]
#        self.comparename  = CheckTools(args).compare[1]
#        self.venn_cluster = CheckTools(args).venn_cluster[0]
#        self.mdspe        = args.get('mdspe')
#        self.keggspe      = args.get('kegg')
        self.BioScripts   = BioConfig.config.get('srnaenv','BioScripts')
        self.python276    = BioConfig.config.get('srnaenv','python_v276')
        self.python2710   = BioConfig.config.get('srnaenv','python_v2710')
        self.blast2226    = BioConfig.config.get('software','blast1')
        self.infernal     = BioConfig.config.get('software','infernal')
#        self.Rfam         = BioConfig.config.get('database','Rfam')
        self.miRBase21    = BioConfig.config.get('database','miRBase21')
#        self.kobas_v1     = BioConfig.config.get('database','kobas_v20120208')
#        self.kobas_v2     = BioConfig.config.get('database','kobas_v20140801')
#        self.kobas        = BioConfig.config.get('database','kobas')
        self.R_v2153      = BioConfig.config.get('srnaenv','R_v2153')
        self.LIBRARY_PATH = BioConfig.config.get('srnaenv','LIBRARY_PATH')
        self.PyPakage     = BioConfig.config.get('srnaenv','Site_Packages')
        self.perlExec     = BioConfig.config.get('srnaenv','perl_v5182')
        self.Sshow        = BioConfig.config.get('srnaenv','Sshow')
        self.miranda      = BioConfig.config.get('software','miranda')
        self.RNAhybrid    = BioConfig.config.get('software','RNAhybrid')
        self.sched        = args.get('sched')
        self.prefix       = args.get('abbr').strip().split(',')[0]            #NAT 预测的前缀名 和靶基因预测文件的前缀
        self.org          = args.get('org')
        self.common       = args.get('common')
        self.project      = args.get('project')
#        self.contract     = args.get('contract')
#        self.gff3         = args.get('gff3')
#        self.length       = args.get('length')
#        self.koann        = args.get('ko')
        self.type         = args.get('type')
#        self.refer        = CheckTools(self.args).getRefConfig().get('refer')
        self.gtf          = CheckTools(self.args).getRefConfig().get('gtf')
#        self.gene         = CheckTools(self.args).getRefConfig().get('gene')
        self.geneAnn      = CheckTools(self.args).getRefConfig().get('geneAnn')
#        self.go           = CheckTools(self.args).getRefConfig().get('go')
        self.utr3         = CheckTools(self.args).getRefConfig().get('utr3')


    def run(self):
#        self.edit_family_analy()

        if self.org.lower() == 'refplant':
            self.mirna_target_refplant()
            if len(self.sample_list) != 1:
                self.kegg_blast_refplant()
                for each_compare in self.comparename.replace(':','vs').split(','):
                    self.GO_and_kegg_enrich_refplant(each_compare)
            else:
                self.GO_and_kegg_enrich_refplant_single(self.sample_list[0]) 

        elif self.org.lower() == 'norefplant':
            self.mirna_target_norefplant()
            if len(self.sample_list) != 1:
                for each_compare in self.comparename.replace(':','vs').split(','):
                    self.GO_and_kegg_enrich_norefplant(each_compare)
            else:
                self.GO_and_kegg_enrich_norefplant_single(self.sample_list[0])

        elif self.org.lower() == 'refanimal':                
            self.mirna_target_refanimal()
#            if len(self.sample_list) != 1:
#                self.kegg_blast_refanimal()
#                for each_compare in self.comparename.replace(':','vs').split(','):
#                    self.GO_and_kegg_enrich_refanimal(each_compare)
#            else:
#                self.GO_and_kegg_enrich_refanimal_single(self.sample_list[0])

        elif self.org.lower() == 'norefanimal':
            self.mirna_target_norefanimal()
            if len(self.sample_list) != 1:
                for each_compare in self.comparename.replace(':','vs').split(','):
                    self.GO_and_kegg_enrich_norefanimal(each_compare)
            else:
                self.GO_and_kegg_enrich_norefanimal_single(self.sample_list[0])

#        if len(self.sample_list) == 1:
#            self.single_diff_analy()
#        else:
#            if self.venn_cluster:
#                self.diff_analy_venn()
#            else:
#                self.diff_analy_no_venn()
#        self.srna_result()
#        self.srna_report()


    def edit_family_analy(self):
       
        edit_analy   =  os.path.join(self.projpath,'11.edit_family','edit_analy')
        family_analy =  os.path.join(self.projpath,'11.edit_family','family_analy')
        make_dir(edit_analy)
        make_dir(family_analy)

        code = (
                'date +"%D %T ->Start edit and family analysis " && \\\n\n'
                'cat {projpath}/3.known/known_miRNAs/hairpin.fa \\\n'
                '    {projpath}/8.novel/novel_miRNAs/hairpin.fa \\\n'
                '    > {edit_analy}/known_novel_hairpin.fa && \n'
                'cat {projpath}/3.known/known_miRNAs/hairpin_mature.fa \\\n'
                '    {projpath}/8.novel/novel_miRNAs/hairpin_mature.fa \\\n'
                '    > {edit_analy}/known_novel_hairpin_mature.fa  && \n\n'
                'export PATH={mirdeep2}:$PATH  && \n'
                'cd {edit_analy} && \n'
        ).format( projpath = self.projpath,mirdeep2 = self.mirdeep2,edit_analy = edit_analy)

        for sampid in self.sample_list:
            code += (
                     '{perlExec} {BioScripts}/10.0edit_family/quantifier_gb.pl \\\n'
                     '      -p {edit_analy}/known_novel_hairpin.fa \\\n'
                     '      -m {edit_analy}/known_novel_hairpin_mature.fa \\\n'
                     '      -r {projpath}/2.map/{sampid}/{sampid}.collapse.fa \\\n'
                     '      -y {sampid} \\\n'
                     '      -g 2 -T 10 -d && \n\n'
                     '{python276} {BioScripts}/10.0edit_family/miRNA_editing.py \\\n'
                     '      {edit_analy}/expression_analyses/expression_analyses_{sampid}/miRBase.mrd \\\n'
                     '      {edit_analy}/{sampid} && \n\n'
                     'awk -F"\\t" \'{{if(/^>/){{if($2){{marker=1}}else{{marker=0}}}}if(marker){{print}}}}\' \\\n' 
                     '      {edit_analy}/{sampid}_editing_stats.txt | head -10 \\\n' 
                     '      > {edit_analy}/{sampid}_editing_stats.example.txt && \n\n'
            ).format(projpath = self.projpath, BioScripts = self.BioScripts,perlExec = self.perlExec,
                     sampid = sampid, python276 = self.python276,edit_analy = edit_analy)

        code += (
                 'export PATH={blast2226}:{infernal}/bin:$PATH && \n'
                 '{perlExec} {BioScripts}/10.0edit_family/rfam_scan.pl \\\n'
                 '      -blastdb {Rfam}/Rfam.fasta \\\n'
                 '           {Rfam}/Rfam.cm \\\n'
                 '           {projpath}/8.novel/novel_miRNAs/hairpin.fa \\\n'
                 '      -bt 1e-05 \\\n'
                 '      -o {family_analy}/novel_hairpin.rfam.gff3 && \n\n'
                 '{python2710} {BioScripts}/10.0edit_family/miRNA_family.py \\\n'
                 '      --known   {projpath}/3.known/known_miRNAs/hairpin.fa \\\n'
                 '      --unknown {projpath}/8.novel/novel_miRNAs/hairpin.fa \\\n'
                 '      --prefix  {prefix} \\\n'
                 '      --family_hp {miRBase21}/hairpin.fa \\\n'
                 '      --family_miFam {miRBase21}/miFam.dat \\\n'
                 '      --novo_gff {family_analy}/novel_hairpin.rfam.gff3 && \n\n'
                 '{perlExec}  {BioScripts}/10.0edit_family/miRNA_family_show.pl \\\n'
                 '      {edit_analy}/{prefix}_miRNA_family.txt \\\n'
                 '      {family_analy}/{prefix}_miRNA_family && \n\n'
                 'cut -f 1-10 {family_analy}/{prefix}_miRNA_family.mir_num.txt | \\\n'
                 '      awk -F"\\t" \'{{if($2>0 || $3>0){{print}}}}\' | head \\\n'
                 '      > {family_analy}/{prefix}_miRNA_family.mir_num.example.txt && \n'
                 'cut -f 1-10 {family_analy}/{prefix}_miRNA_family.mir_sign.txt | \\\n'
                 '      awk -F"\\t" \'{{if($2=="+" || $3=="+"){{print}}else if(NR==1){{print}}}}\' | head \\\n'
                 '      > {family_analy}/{prefix}_miRNA_family.mir_sign.example.txt && \n\n'
                 'date +"%D %T ->Finish edit and family analysis "\n'
        ).format(blast2226 = self.blast2226, infernal = self.infernal, perlExec = self.perlExec,
                 BioScripts = self.BioScripts, Rfam = self.Rfam, python2710 = self.python2710,
                 projpath = self.projpath, prefix = self.prefix, miRBase21 = self.miRBase21,
                 edit_analy = edit_analy, family_analy = family_analy)

        shell_path = os.path.join(self.projpath,'11.edit_family','edit_family_analy.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'edit_family_analy', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='novel_miRNA_stat')
        self.job.add_order(job_name, after_jobs='srna_result')
   
    def mirna_target_refplant(self):
        
        code = (
                'set -e \\\n'
                'date +"%D %T ->Start refplant mirna target predicted " && \\\n'
                'cd {projpath}/12.target && \\\n'
                'ln -sf {gene} \\\n'
                '       {projpath}/12.target && \\\n'
                'sed \'s/u/t/g\' {projpath}/8.novel/novel_miRNAs/mature.fa \\\n'
                '    > {projpath}/12.target/novel_mature.fa && \n'
                'cat {projpath}/3.known/known_miRNAs/mature.fa \\\n'
                '    {projpath}/12.target/novel_mature.fa \\\n'
                '    > {projpath}/12.target/mature.fa && \\\n'
                'rm -rf {projpath}/12.target/novel_mature.fa && \n\n'
                'echo "植物靶基因预测软件更换为在线预测软件psRNATarget,需要将{projpath}/12.target/mature.fa" \n'
                'echo "和{gene}下载下来" \n'
                'echo "在http://plantgrn.noble.org/psRNATarget/analysis?function=3网站上输入文件进行在线预测" \n'
                'echo "将结果文件下载下来后上传到{projpath}/12.target路径下重新投递job即可" \n\n'
                'mv {projpath}/12.target/psRNATargetJob*.txt  \\\n'
                '   {projpath}/12.target/{abbr}_targets.txt && \n'
                'ln -sf {abbr}_targets.txt {projpath}/12.target/{abbr}_targets && \n'
                'awk \'BEGIN{{FS=OFS="\\t"}}{{print $1,$2}}\' {projpath}/12.target/{abbr}_targets.txt \\\n'
                '     > {projpath}/12.target/{abbr}_targets_trans.txt && \n'
                '{perlExec}  {BioScripts}/10.2target/trans_gene_pairs.pl \\\n'
                '     -gtf {gtf} \\\n'
                '     -tp {projpath}/12.target/{abbr}_targets \\\n'
                '     -o  {projpath}/12.target/{abbr}_targets.pairs && \n'
                'awk \'{{if(name!=$1){{print}}}}\' {projpath}/12.target/{abbr}_targets.pairs | \\\n'
                '     awk \'{{print $1"\\t"$2"\\t"$3}}\' | head -15 \\\n'
                '     > {projpath}/12.target/{abbr}_targets.pairs_example.xls && \n'
                'awk \'{{print $1"\\t"$3}}\' {projpath}/12.target/{abbr}_targets.pairs | sort -u \\\n'
                '     > {projpath}/12.target/{abbr}_targets_gene.pairs && \n'
                '{perlExec}  {BioScripts}/10.2target/targets.pairs_transfer_v4.pl \\\n'
                '     {projpath}/12.target/{abbr}_targets_gene.pairs \\\n'
                '     {geneAnn} \\\n'
                '     {projpath}/12.target/{abbr}_targets.pairs && \\\n\n'
                'date +"%D %T ->Finish refplant mirna target predicted "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 gene = self.gene,abbr = self.prefix,gtf = self.gtf,geneAnn = self.geneAnn)

        shell_path = os.path.join(self.projpath,'12.target','mirna_target_refplant.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'mirna_target_refplant', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='novel_miRNA_stat')


    def mirna_target_norefplant(self):
        
        code = (
                'set -e \\\n'
                'date +"%D %T ->Start refplant mirna target predicted " && \\\n'
                'cd {projpath}/12.target && \\\n'
                'ln -sf {gene} \\\n'
                '       {projpath}/12.target && \\\n'
                'sed \'s/u/t/g\' {projpath}/8.novel/novel_miRNAs/mature.fa \\\n'
                '    > {projpath}/12.target/novel_mature.fa && \n'
                'cat {projpath}/3.known/known_miRNAs/mature.fa \\\n'
                '    {projpath}/12.target/novel_mature.fa \\\n'
                '    > {projpath}/12.target/mature.fa && \\\n'
                'rm -rf {projpath}/12.target/novel_mature.fa && \n\n'
                'echo "植物靶基因预测软件更换为在线预测软件psRNATarget,需要将{projpath}/12.target/mature.fa" \n'
                'echo "和{gene}下载下来" \n'
                'echo "在http://plantgrn.noble.org/psRNATarget/analysis?function=3网站上输入文件进行在线预测" \n'
                'echo "将结果文件下载下来后上传到{projpath}/12.target路径下重新投递job即可" \n\n'
                'mv {projpath}/12.target/psRNATargetJob*.txt  \\\n'
                '   {projpath}/12.target/{abbr}_targets.txt && \n'
                'ln -sf {abbr}_targets.txt {projpath}/12.target/{abbr}_targets && \n'
                'awk \'BEGIN{{FS=OFS="\\t"}}{{print $1,$2}}\' {projpath}/12.target/{abbr}_targets.txt \\\n'
                '    > {projpath}/12.target/{abbr}_targets_trans.txt && \n'
                'sed -i /^#.*/d {projpath}/12.target/{abbr}_targets_trans.txt && \n'
                'cp {projpath}/12.target/{abbr}_targets_trans.txt \\\n'
                '   {projpath}/12.target/{abbr}_targets.pairs && \n'
                'awk \'BEGIN{{print "miRNA\\ttarget_mRNA"}}{{if(name!=$1){{print}}}}\' \\\n'
                '    {projpath}/12.target/{abbr}_targets.pairs | head -15 \\\n'
                '    >{projpath}/12.target/{abbr}_targets_example.xls && \n'
                'awk \'{{print $1"\\t"$2}}\' {projpath}/12.target/{abbr}_targets.pairs | sort -u \\\n'
                '    >{projpath}/12.target/{abbr}_targets_gene.pairs && \n'
                '{perlExec}  {BioScripts}/10.2target/targets.pairs_transfer_v2.pl \\\n'
                '     {projpath}/12.target/{abbr}_targets_gene.pairs \\\n'
                '     {geneAnn} \\\n\n'
                'date +"%D %T ->Finish refplant mirna target predicted "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 gene = self.gene,abbr = self.prefix,gtf = self.gtf,geneAnn = self.geneAnn)

        shell_path = os.path.join(self.projpath,'12.target','mirna_target_norefplant.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'mirna_target_norefplant', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='novel_miRNA_stat')

    def mirna_target_refanimal(self):

        code = (
                'date +"%D %T ->Start mirna target split " && \\\n'
                'cd {projpath}/12.target && \\\n'
                'ln -sf {utr3} {projpath}/12.target && \\\n'
                'cat {projpath}/3.known/known_miRNAs/mature.fa \\\n'
#                '    {projpath}/8.novel/novel_miRNAs/mature.fa \\\n'
                '    > {projpath}/12.target/mature.fa && \\\n'
                '{pythonExec} {BioScripts}/10.2target/split_fa_by_fixnum.py \\\n'
                '    {projpath}/12.target/mature.fa 40 && \\\n'
                'date +"%D %T ->Start mirna target split "'
                ).format(projpath = self.projpath,BioScripts = self.BioScripts,
                         pythonExec = self.python2710,utr3 = self.utr3)
         
        shell_path = os.path.join(self.projpath,'12.target','mirna_target_split.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'mirna_target_split', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='known_miRNA_stat')

        for num in xrange(1,5):
            code_miRanda = (
                            'date +"%D %T ->Start mirna target predicted by miranda " && \\\n\n'
                            '{miranda}  {projpath}/12.target/mature.fa_{num}.fasta  \\\n'
                            '   {utr3} \\\n'
                            '   -sc 140 -en -10 -scale 4 -strict \\\n'
                            '   -out  {projpath}/12.target/miRanda/miranda_targets_out_{num} && \n\n' 
                            '{perlExec} {BioScripts}/10.2target/mi_result_fmt.pl \\\n'
                            '   {projpath}/12.target/miRanda/miranda_targets_out_{num} \\\n'
                            '   {projpath}/12.target/miRanda/miranda_targets_out_{num}.fmt && \n\n' 
                            '{perlExec}  {BioScripts}/10.2target/mi_result.pl \\\n'
                            '   -i {projpath}/12.target/miRanda/miranda_targets_out_{num} \\\n'
                            '   -o {projpath}/12.target/miRanda/miranda_targets_{num} && \n\n'
                            'awk \'{{print substr($1,3)"\\t"$2}}\' {projpath}/12.target/miRanda/miranda_targets_{num} | '
                            'sort -u > {projpath}/12.target/miRanda/miranda_targets_{num}.pairs && \n\n'
                            'date +"%D %T ->Finish mirna target predicted by miranda " \n'
                           ).format(miranda = self.miranda,projpath = self.projpath, num = num,
                                    utr3 = self.utr3,perlExec = self.perlExec,BioScripts = self.BioScripts)
            shell_path = os.path.join(self.projpath,'12.target','miRanda','miranda_runtarget_{num}.sh'.format(num = num))
            write_shell(code_miRanda,shell_path)
            job_name = os.path.splitext(os.path.basename(shell_path))[0]
            atom_job = AtomJob(job_name, shell_path,1,'miranda_runtarget', sched=self.sched)
            self.job.add_atom_job(atom_job)
            self.job.add_order(job_name, previous_jobs='mirna_target_split')
            self.job.add_order(job_name, after_jobs='target_Common')

            code_pita = (
                         'date +"%D %T ->Start mirna target predicted by PITA " && \\\n'
                         'cd  {projpath}/12.target/PITA && \\\n' 
                         '{perlExec}  {BioScripts}/10.2target/pita_prediction.pl \\\n'
                         '  -utr {utr3} \\\n'
                         '  -mir {projpath}/12.target/mature.fa_{num}.fasta \\\n'
                         '  -prefix PITA_{num} && \\\n'
                         'date +"%D %T ->Finish mirna target predicted by PITA "\n'
                        ).format(projpath = self.projpath,perlExec = self.perlExec,BioScripts = self.BioScripts,
                                 utr3 = self.utr3, num = num)

            shell_path = os.path.join(self.projpath,'12.target','PITA','pita_runtarget_{num}.sh'.format(num = num))
            write_shell(code_pita,shell_path)
            job_name = os.path.splitext(os.path.basename(shell_path))[0]
            atom_job = AtomJob(job_name, shell_path,1,'pita_runtarget', sched=self.sched)
            self.job.add_atom_job(atom_job)
            self.job.add_order(job_name, previous_jobs='mirna_target_split')
            self.job.add_order(job_name, after_jobs='target_Common')

            code_RNAhybrid = (
                              'date +"%D %T ->Start mirna target predicted by RNAhybrid " && \\\n'
                              '{RNAhybrid} -s {type} \\\n'
                              ' -t {utr3} \\\n'
                              ' -q {projpath}/12.target/mature.fa_{num}.fasta \\\n'
                              ' -e -10 -p 0.05 \\\n'
                              ' > {projpath}/12.target/RNAhybrid/RNAhybrid_miRNA_target_{num}_pairs && \\\n'
                              'date +"%D %T ->Finish mirna target predicted by RNAhybrid "\n'
                             ).format(projpath = self.projpath,RNAhybrid = self.RNAhybrid,type = self.type,
                                      utr3 = self.utr3,num = num)

            shell_path = os.path.join(self.projpath,'12.target','RNAhybrid','RNAhybrid_runtarget_{num}.sh'.format(num = num))
            write_shell(code_RNAhybrid,shell_path)
            job_name = os.path.splitext(os.path.basename(shell_path))[0]
            atom_job = AtomJob(job_name, shell_path,1,'RNAhybrid_runtarget', sched=self.sched)
            self.job.add_atom_job(atom_job)
            self.job.add_order(job_name, previous_jobs='mirna_target_split')
            self.job.add_order(job_name, after_jobs='target_Common')

        
        code_Common = (
                       'date +"%D %T ->Start target information integration " && \\\n'
                       'cat {projpath}/12.target/miRanda/miranda_targets_*.pairs \\\n'
                       '    > {projpath}/12.target/miRanda/miranda_targets.pairs && \n'
                       'cat {projpath}/12.target/miRanda/miranda_targets_out_*.fmt \\\n'
                       '    > {projpath}/12.target/miRanda/miranda_targets_out.fmt && \n'
                       'cat {projpath}/12.target/PITA/PITA_*_pita_results_targets.tab \\\n'
                       '    > {projpath}/12.target/PITA/PITA_pita_results_targets.tab && \n'
                       'cat {projpath}/12.target/RNAhybrid/RNAhybrid_miRNA_target_*_pairs \\\n'
                       '    > {projpath}/12.target/RNAhybrid/RNAhybrid_miRNA_target_pairs && \n\n'
                       'cd {projpath}/12.target/Common && \n'
                       '{perlExec}  {BioScripts}/10.2target/get_RNAhybrid_PITA_miRanda.pl \\\n'
                       '    -RNAhybrid {projpath}/12.target/RNAhybrid/RNAhybrid_miRNA_target_pairs \\\n'
                       '    -PITA      {projpath}/12.target/PITA/PITA_pita_results_targets.tab \\\n'
                       '    -miRanda   {projpath}/12.target/miRanda/miranda_targets.pairs && \n\n'
                       '{perlExec}  {BioScripts}/10.2target/trans_gene_pairs_rad.pl \\\n'
                       '    -gtf {gtf} \\\n'
                       '    -tp  {projpath}/12.target/Common/all_target.xls \\\n'
                       '    -o   {projpath}/12.target/Common/all_target_gene.xls && \n\n'
                       'head -15 {projpath}/12.target/Common/commom_target.xls \\\n'
                       '    > {projpath}/12.target/Common/commom_target_example.xls && \n'
                       'awk \'{{print $1"\\t"$3}}\' {projpath}/12.target/Common/all_target_gene.xls | \\\n'
                       '    sort -u > {projpath}/12.target/Common/all_targets_gene.pairs && \n'
                       'ln -sf {projpath}/12.target/Common/all_targets_gene.pairs \\\n'
                       '    {projpath}/12.target/{abbr}_targets_gene.pairs && \n'
                       '{perlExec}  {BioScripts}/10.2target/targets.pairs_transfer_v4.pl \\\n'
                       '    {projpath}/12.target/Common/all_targets_gene.pairs \\\n'
                       '    {geneAnn} \\\n'
                       '    {projpath}/12.target/Common/all_target_gene.xls && \n'
                       '{pythonExec} {BioScripts}/10.2target/target_clean.py \\\n'
                       '    --outdir {projpath}/12.target && \n'
                       'date +"%D %T ->Finish target information integration "\n'
                      ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                               pythonExec = self.python2710,geneAnn = self.geneAnn,gtf = self.gtf,abbr = self.prefix) 
           
        shell_path = os.path.join(self.projpath,'12.target','Common','target_Common.sh')
        write_shell(code_Common,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'target_Common', sched=self.sched)
        self.job.add_atom_job(atom_job)

    def mirna_target_norefanimal(self):
        
        code = (
                'date +"%D %T ->Start mirna target norefanimal " && \\\n'
                'cd {projpath}/12.target && \n '
                'ln -sf {refer} {projpath}/12.target && \n'
                'ln -sf {gff3}  {projpath}/12.target && \n'
                'cat {projpath}/3.known/known_miRNAs/mature.fa \\\n'
                '    {projpath}/8.novel/novel_miRNAs/mature.fa \\\n'
                '    > {projpath}/12.target/mature.fa && \\\n'
                '{perlExec} {BioScripts}/10.2target/Trinity_gff3-3UTR.pl \\\n'
                '   {refer} \\\n'
                '   {gff3} 16 \\\n'
                '   > {projpath}/12.target/Trinity_3UTR.fasta && \n'
                '{miranda}  {projpath}/12.target/mature.fa \\\n'
                '   {projpath}/12.target/Trinity_3UTR.fasta \\\n'
                '   -sc 140 -en -10 -scale 4 -strict \\\n'
                '   -out {projpath}/12.target/mature.miranda_targets_out && \n'
                '{perlExec} {BioScripts}/10.2target//mi_result_fmt.pl \\\n'
                '   {projpath}/12.target/mature.miranda_targets_out \\\n'
                '   {projpath}/12.target/mature.miranda_targets_out.fmt && \n'
                '{perlExec} {BioScripts}/10.2target/mi_result.pl \\\n'
                '   -i {projpath}/12.target/mature.miranda_targets_out \\\n'
                '   -o {projpath}/12.target/mature.miranda_targets && \n'
                'awk \'BEGIN{{print "miRNA\\ttarget_gene"}}{{if(name!=$1){{print}}}}\' \\\n'
                '   {projpath}/12.target/mature.miranda_targets.pairs | \\\n'
                '   head -15 > {projpath}/12.target/mature.miranda_targets.pairs.example && \n'
                'ln -sf {projpath}/12.target/mature.miranda_targets.pairs.example \\\n'
                '   {projpath}/12.target/mature.miranda_targets_example.xls && \n'
                'ln -sf {projpath}/12.target/mature.miranda_targets.pairs \\\n'
                '   {projpath}/12.target/mature.miranda_targets_gene.pairs && \n'
                '{perlExec} {BioScripts}/10.2target/targets.pairs_transfer_v2.pl \\\n'
                '   {projpath}/12.target/mature.miranda_targets.pairs \\\n'
                '   {geneAnn} && \n'
                'date +"%D %T ->Finish mirna target norefanimal "\n'
               ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                        refer = self.refer,gff3 = self.gff3,miranda = self.miranda,geneAnn = self.geneAnn)
 
        shell_path = os.path.join(self.projpath,'12.target','mirna_target_norefanimal.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'mirna_target_norefanimal', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='novel_miRNA_stat')


    def diff_analy_venn(self):
        
        diffdir =  os.path.join(self.projpath,'13.diff','diffAnalysisResult')
        make_dir(diffdir)
     
        code = (
                'cd {projpath}/13.diff && \\\n'
                '{perlExec}  {BioScripts}/10.1diff_aly/run_DE_sRNA_v2.pl \\\n'
                '   -r  {projpath}/8.novel/novel_miRNAs/mature.readcount,{projpath}/3.known/known_miRNAs/mature.readcount \\\n'
                '   -s          {sample} \\\n'
                '   -group      {group} \\\n'
                '   -groupname  {groupname} \\\n'
                '   -g          {compare} \\\n'
                '   -venn       {venn_cluster} \\\n'
                '   -o          {diffdir} '
                '   > {projpath}/13.diff/diffAnalysis.sh '
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 diffdir = diffdir,sample = self.sample,group = self.group,groupname = self.groupname,
                 compare = self.compare, venn_cluster = self.venn_cluster)
        
        assert not os.system(code)
        shell_path = os.path.join(self.projpath,'13.diff','diffAnalysis.sh')
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'diffAnalysis', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs = 'novel_miRNA_stat')

    def diff_analy_no_venn(self):
    
        diffdir =  os.path.join(self.projpath,'13.diff','diffAnalysisResult')
        make_dir(diffdir)
     
        code = (
                'cd {projpath}/13.diff && \\\n'
                '{perlExec}  {BioScripts}/10.1diff_aly/run_DE_sRNA_v2.pl \\\n'
                '   -r  {projpath}/8.novel/novel_miRNAs/mature.readcount,{projpath}/3.known/known_miRNAs/mature.readcount \\\n'
                '   -s          {sample} \\\n'
                '   -group      {group} \\\n'
                '   -groupname  {groupname} \\\n'
                '   -g          {compare} \\\n'
                '   -o          {diffdir} '
                '   > {projpath}/13.diff/diffAnalysis.sh '
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 diffdir = diffdir,sample = self.sample,group = self.group,groupname = self.groupname,
                 compare = self.compare)
        
        assert not os.system(code)
        shell_path = os.path.join(self.projpath,'13.diff','diffAnalysis.sh')
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'diffAnalysis', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs = 'novel_miRNA_stat')
    

    def single_diff_analy(self):
                
        diffdir =  os.path.join(self.projpath,'13.diff','diffAnalysisResult')
        make_dir(diffdir)
     
        code = (
                'cd {projpath}/13.diff && \\\n'
                '{perlExec}  {BioScripts}/10.1diff_aly/run_DE_single_sRNA_v2.pl \\\n'
                '   -r  {projpath}/8.novel/novel_miRNAs/mature.readcount,{projpath}/3.known/known_miRNAs/mature.readcount \\\n'
                '   -s          {sample} \\\n'
                '   -group      {group} \\\n'
                '   -groupname  {groupname} \\\n'
                '   -g          {compare} \\\n'
                '   -o          {diffdir} '
                '   > {projpath}/13.diff/diffAnalysis.sh '
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 diffdir = diffdir,sample = self.sample,group = self.group,groupname = self.groupname,
                 compare = self.compare)
        
        assert not os.system(code)
        shell_path = os.path.join(self.projpath,'13.diff','diffAnalysis.sh')
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'diffAnalysis', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs = 'novel_miRNA_stat')


    def kegg_blast_refplant(self):
        '''
           运行此步的条件:有参植物多样本
        '''

        code = (
                'date +"%D %T ->Start refplant kegg blast " && \\\n'
                'cd {projpath}/Blast && \\\n\n'
                'PYTHONPATH={kobas_v20140801}/src:{kobas}/src:$PYTHONPATH \\\n'
                'sort -k 1,1 {projpath}/13.diff/ALL.DElist.txt \\\n'
                '   > {projpath}/Blast/diffmiRNAID && \n\n'
                'sort -k 1,1 {projpath}/12.target/{abbr}_targets_gene.pairs | \\\n'
                '   join {projpath}/Blast/diffmiRNAID - | sort -u \\\n'
                '   > {projpath}/Blast/diffmiRNA-gene.pairs && \n'
                'awk \'{{print $2}}\' {projpath}/Blast/diffmiRNA-gene.pairs | sort -u \\\n'
                '   > {projpath}/Blast/diffmiRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/runKEGG.pl \\\n'
                '   -directory-gtf {gtf} \\\n'
                '   -gene-id {projpath}/Blast/diffmiRNA-geneid \\\n'
                '   -exist-genome y \\\n'
                '   -directory-kobas-output {projpath}/Blast \\\n'
                '   -directory-genome {refer} \\\n'
                '   -species {keggspe} \\\n'
                '   -samplename KEGG  && \\\n\n'
                'date +"%D %T ->Finish refplant kegg blast "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 abbr = self.prefix,refer = self.refer,gtf = self.gtf,kobas_v20140801 = self.kobas_v2,
                 kobas = self.kobas,keggspe = self.keggspe)

        shell_path = os.path.join(self.projpath,'Blast','kegg_blast_refplant.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,4,'kegg_blast_refplant', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_refplant','diffAnalysis'])


    def kegg_blast_refanimal(self):
        '''
            运行此步的条件:有参动物多样本
        '''
        
        code = (
                'date +"%D %T ->Start refanimal kegg blast " && \\\n'
                'cd {projpath}/Blast && \\\n\n'
                'PYTHONPATH={kobas_v20140801}/src:{kobas}/src:$PYTHONPATH \\\n'
                'sort -k 1,1 {projpath}/13.diff/ALL.DElist.txt \\\n'
                '   > {projpath}/Blast/diffmiRNAID && \n\n'
                'sort -k 1,1 {projpath}/12.target/Common/all_targets_gene.pairs | \\\n'
                '   join {projpath}/Blast/diffmiRNAID - | sort -u \\\n'
                '   > {projpath}/Blast/diffmiRNA-gene.pairs && \n'
                'awk \'{{print $2}}\' {projpath}/Blast/diffmiRNA-gene.pairs | sort -u \\\n'
                '   > {projpath}/Blast/diffmiRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/runKEGG.pl \\\n'
                '   -directory-gtf {gtf} \\\n'
                '   -gene-id {projpath}/Blast/diffmiRNA-geneid \\\n'
                '   -exist-genome y \\\n'
                '   -directory-kobas-output {projpath}/Blast \\\n'
                '   -directory-genome {refer} \\\n'
                '   -species {keggspe} \\\n'
                '   -samplename KEGG  && \\\n\n'
                'date +"%D %T ->Finish refplant kegg blast "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 abbr = self.prefix,refer = self.refer,gtf = self.gtf,kobas_v20140801 = self.kobas_v2,
                 kobas = self.kobas,keggspe = self.keggspe)

        shell_path = os.path.join(self.projpath,'Blast','kegg_blast_refanimal.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,4,'kegg_blast_animal', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['target_Common','diffAnalysis'])


    def GO_and_kegg_enrich_refplant(self,each_compare):
                 
        code = (
                'date +"%D %T ->Start GO and kegg enrichment " && \\\n\n'
                'cd {projpath}/14.enrich/{each_compare} && \\\n'
                'sort -k 1,1 {projpath}/13.diff/diffAnalysisResult/{each_compare}/{each_compare}.DElist.txt \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID && \n\n'
                'sort -k 1,1 {projpath}/12.target/{abbr}_targets_gene.pairs | \\\n'
                '   join {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID - | sort -u \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs && \n\n'
                '{perlExec} {BioScripts}/11Enrichment/add_description.pl \\\n'
                '   {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs \\\n'
                '   {geneAnn} && \n\n'
                'mv {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs.mg \\\n'
                '   {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs && \n\n'
                'awk \'{{if(NR>0){{print $2}}}}\' {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs |sort -u \\\n'
                '   >{projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph_v3.pl \\\n'
                '   -i {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -goann {go} \\\n'
                '   -gtf {gtf} \\\n'
                '   -o {projpath}/14.enrich/{each_compare}/GO2 \\\n'
                '   -p {each_compare} && \n\n'
                'mkdir -p {projpath}/14.enrich/{each_compare}/Pathway && \n'
                'export PATH={R_v2153}:$PATH && \n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/runKEGG_v1.pl \\\n'
                '   -directory-gtf {gtf} \\\n'
                '   -gene-id {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -exist-genome y \\\n'
                '   -directory-kobas-output {projpath}/14.enrich/{each_compare}/Pathway \\\n'
                '   -directory-genome {refer} \\\n'
                '   -species {keggspe} \\\n'
                '   -samplename {each_compare} && \n\n'
                'cd {projpath}/14.enrich/{each_compare}/Pathway && \n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table add.{each_compare}.identify.xls \\\n' 
                '   --abbr {keggspe} && \\\n\n'
                'date +"%D %T ->Finish GO and kegg enrichment "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,each_compare = each_compare,gtf = self.gtf,geneAnn = self.geneAnn,
                 abbr = self.prefix,refer = self.refer,R_v2153 = self.R_v2153,kobas_v1 = self.kobas_v1,
                 kobas_v2 = self.kobas_v2, LIBRARY_PATH = self.LIBRARY_PATH, kobas = self.kobas,
                 PyPakage = self.PyPakage, pythonExec = self.python2710,keggspe = self.keggspe)
         
        shell_path = os.path.join(self.projpath,'14.enrich',each_compare,'{}_GO_and_kegg_enrich.sh'.format(each_compare))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='kegg_blast_refplant')
        self.job.add_order(job_name, after_jobs='srna_result')


    def GO_and_kegg_enrich_refplant_single(self,sample):
        
        code = (
                'date +"%D %T ->Start single sample GO and kegg enrich analysis " && \\\n\n'
                'cd {projpath}/14.enrich && \\\n'
                'awk \'{{print $2}}\' {projpath}/12.target/{abbr}_targets_gene.pairs | sort -u \\\n' 
                '     >{projpath}/14.enrich/{abbr}.miRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph.pl \\\n'
                '   -i     {projpath}/14.enrich/{abbr}.miRNA-geneid \\\n'
                '   -goann {go} \\\n'
                '   -gtf   {gtf} \\\n'
                '   -o     {projpath}/14.enrich/GO_enrichment \\\n'
                '   -p     {sample} && \n\n'
                'sort -k 1,1 {go} | join {projpath}/14.enrich/{abbr}.miRNA-geneid - | sed \'s/ /\\t/g\''
                '>{projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt && \\\n'
                'mkdir -p {projpath}/14.enrich/GO_classify && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/go_classification_v3.pl \\\n'
                '   {projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt \\\n'
                '   {projpath}/14.enrich/GO_classify && \n\n'
                'mkdir -p {projpath}/14.enrich/Pathway && \\\n'
                'export PATH={R_v2153}:$PATH && \\\n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \\\n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/runKEGG.pl \\\n'
                '   -directory-gtf {gtf} \\\n'
                '   -gene-id {projpath}/14.enrich/{abbr}.miRNA-geneid \\\n '
                '   -exist-genome y \\\n'
                '   -directory-kobas-output {projpath}/14.enrich/Pathway \\\n'
                '   -directory-genome {refer} \\\n'
                '   -species {keggspe} \\\n'
                '   -samplename {sample} && \n\n'
                'cd {projpath}/14.enrich/Pathway && \\\n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table {projpath}/14.enrich/Pathway/add.${sample}.identify.xls \\\n'
                '   --abbr {keggspe} && \n\n'
                'date +"%D %T ->Finish "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,gtf = self.gtf,sample = sample,abbr = self.prefix,refer = self.refer,
                 pythonExec = self.python2710, R_v2153 = self.R_v2153,LIBRARY_PATH = self.LIBRARY_PATH,
                 PyPakage = self.PyPakage, kobas_v1 = self.kobas_v1, kobas_v2 = self.kobas_v2, kobas =self.kobas,
                 keggspe = self.keggspe)
        
        shell_path = os.path.join(self.projpath,'14.enrich','{}_GO_and_kegg_enrich_refplant_single.sh'.format(sample))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_refplant','diffAnalysis'])
        self.job.add_order(job_name, after_jobs='srna_result')
        
    def GO_and_kegg_enrich_norefplant(self,each_compare):
        
        code = (
                'date +"%D %T ->Start GO and kegg enrichment " && \\\n\n'
                'cd {projpath}/14.enrich/{each_compare} && \\\n'
                'sort -k 1,1 {projpath}/13.diff/diffAnalysisResult/{each_compare}/{each_compare}.DElist.txt \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID && \n\n'
                'sort -k 1,1 {projpath}/12.target/{abbr}_targets_gene.pairs | \\\n'
                '   join {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID - | sort -u \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs && \n\n'
                'awk \'{{print $2}}\' {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs |sort -u \\\n'
                '   >{projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph_v3.pl \\\n'
                '   -i {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -goann {go} \\\n'
                '   -length {length} \\\n'
                '   -o {projpath}/14.enrich/{each_compare}/GO2 \\\n'
                '   -p {each_compare} && \n\n'
                'mkdir -p {projpath}/14.enrich/GO_classify && \\\n'
                'sort -k 1,1 {go} | \\\n'
                '   join {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid - | \\\n'
                '   sed \'s/ /\\t/g\' > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid.go.txt && \n'
                '{perlExec} {BioScripts}/11Enrichment/go_classification_v3.pl \\\n'
                '   {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid.go.txt \\\n'
                '   {projpath}/14.enrich/each_compare/GO_classify && \n\n'
                'mkdir -p {projpath}/14.enrich/{each_compare}/Pathway && \n'
                'export PATH={R_v2153}:$PATH && \n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/runKEGG_enrich_v1.1.pl  \\\n'
                '   -diff {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -ko   {koann} \\\n'
                '   -g    {each_compare} \\\n'
                '   -t     plant  && \n\n'
                'cd {projpath}/14.enrich/{each_compare}/Pathway && \n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table {projpath}/14.enrich/{each_compare}/{each_compare}.DEG_KEGG_pathway_enrichment_add.xls \\\n' 
                '   --abbr  ko && \\\n\n'
                'date +"%D %T ->Finish GO and kegg enrichment "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,each_compare = each_compare,length = self.length,koann = self.koann,
                 abbr = self.prefix,R_v2153 = self.R_v2153,kobas_v1 = self.kobas_v1,
                 kobas_v2 = self.kobas_v2, LIBRARY_PATH = self.LIBRARY_PATH, kobas = self.kobas,
                 PyPakage = self.PyPakage, pythonExec = self.python2710)

        shell_path = os.path.join(self.projpath,'14.enrich',each_compare,'{}_GO_and_kegg_enrich.sh'.format(each_compare))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_norefplant','diffAnalysis'])
        self.job.add_order(job_name, after_jobs='srna_result')
        

    def GO_and_kegg_enrich_norefplant_single(self,sample):
                
        code = (
                'date +"%D %T ->Start single sample GO and kegg enrich analysis " && \\\n\n'
                'cd {projpath}/14.enrich && \\\n'
                'awk \'{{print $2}}\' {projpath}/12.target/{abbr}_targets_gene.pairs | sort -u \\\n' 
                '     >{projpath}/14.enrich/{abbr}.miRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph.pl \\\n'
                '   -i      {projpath}/14.enrich/{abbr}.miRNA-geneid \\\n'
                '   -goann  {go} \\\n'
                '   -length {length} \\\n'
                '   -o      {projpath}/14.enrich/GO_enrichment \\\n'
                '   -p      {sample} && \n\n'
                'sort -k 1,1 {go} | join {projpath}/14.enrich/{abbr}.miRNA-geneid - | sed \'s/ /\\t/g\''
                '>{projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt && \\\n'
                'mkdir -p {projpath}/14.enrich/GO_classify && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/go_classification_v3.pl \\\n'
                '   {projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt \\\n'
                '   {projpath}/14.enrich/GO_classify && \n\n'
                'mkdir -p {projpath}/14.enrich/Pathway && \\\n'
                'export PATH={R_v2153}:$PATH && \\\n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \\\n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/runKEGG_enrich_v1.1.pl \\\n'               
                '   -diff {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -ko   {koann} \\\n'
                '   -g    {sample} \\\n'
                '   -t     plant  && \n\n'
                'cd {projpath}/14.enrich/Pathway && \\\n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table {projpath}/14.enrich/Pathway/add.${sample}.identify.xls \\\n'
                '   --abbr ko && \n\n'
                'date +"%D %T ->Finish "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,length = self.length,sample = sample,abbr = self.prefix,koann = self.koann,
                 pythonExec = self.python2710, R_v2153 = self.R_v2153,LIBRARY_PATH = self.LIBRARY_PATH,
                 PyPakage = self.PyPakage, kobas_v1 = self.kobas_v1, kobas_v2 = self.kobas_v2, kobas =self.kobas)

        shell_path = os.path.join(self.projpath,'14.enrich','{}_GO_and_kegg_enrich_norefplant_single.sh'.format(sample))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_norefplant','diffAnalysis'])
        self.job.add_order(job_name, after_jobs='srna_result')
        

    def GO_and_kegg_enrich_refanimal(self,each_compare):
        
        code = (
                'date +"%D %T ->Start GO and kegg enrichment " && \\\n\n'
                'cd {projpath}/14.enrich/{each_compare} && \\\n'
                'sort -k 1,1 {projpath}/13.diff/diffAnalysisResult/{each_compare}/{each_compare}.DElist.txt \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID && \n\n'
                'sort -k 1,1 {projpath}/12.target/{abbr}_targets_gene.pairs | \\\n'
                '   join {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID - | sort -u \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs && \n\n'
                '{perlExec} {BioScripts}/11Enrichment/add_description.pl \\\n'
                '   {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs \\\n'
                '   {geneAnn} && \n\n'
                'mv {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs.mg \\\n'
                '   {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs && \n\n'
                'awk \'{{if(NR>0){{print $2}}}}\' {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs |sort -u \\\n'
                '   >{projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph_v3.pl \\\n'
                '   -i {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -goann {go} \\\n'
                '   -gtf {gtf} \\\n'
                '   -o {projpath}/14.enrich/{each_compare}/GO2 \\\n'
                '   -p {each_compare} && \n\n'
                'mkdir -p {projpath}/14.enrich/{each_compare}/Pathway && \n'
                'export PATH={R_v2153}:$PATH && \n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/runKEGG_v1.pl \\\n'
                '   -directory-gtf {gtf} \\\n'
                '   -gene-id {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -exist-genome y \\\n'
                '   -directory-kobas-output {projpath}/14.enrich/{each_compare}/Pathway \\\n'
                '   -directory-genome {refer} \\\n'
                '   -species {keggspe} \\\n'
                '   -samplename {each_compare} && \n\n'
                'cd {projpath}/14.enrich/{each_compare}/Pathway && \n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table add.{each_compare}.identify.xls \\\n' 
                '   --abbr {keggspe} && \\\n\n'
                'date +"%D %T ->Finish GO and kegg enrichment "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,each_compare = each_compare,gtf = self.gtf,geneAnn = self.geneAnn,
                 abbr = self.prefix,refer = self.refer,R_v2153 = self.R_v2153,kobas_v1 = self.kobas_v1,
                 kobas_v2 = self.kobas_v2, LIBRARY_PATH = self.LIBRARY_PATH, kobas = self.kobas,
                 PyPakage = self.PyPakage, pythonExec = self.python2710,keggspe = self.keggspe)
         
        shell_path = os.path.join(self.projpath,'14.enrich',each_compare,'{}_GO_and_kegg_enrich.sh'.format(each_compare))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='kegg_blast_refanimal')
        self.job.add_order(job_name, after_jobs='srna_result')
                

    def GO_and_kegg_enrich_refanimal_single(self,sample):
        
        code = (
                'date +"%D %T ->Start single sample GO and kegg enrich analysis " && \\\n\n'
                'cd {projpath}/14.enrich && \\\n'
                'awk \'{{print $2}}\' {projpath}/12.target/{abbr}_targets_gene.pairs | sort -u \\\n' 
                '     >{projpath}/14.enrich/{abbr}.miRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph.pl \\\n'
                '   -i     {projpath}/14.enrich/{abbr}.miRNA-geneid \\\n'
                '   -goann {go} \\\n'
                '   -gtf   {gtf} \\\n'
                '   -o     {projpath}/14.enrich/GO_enrichment \\\n'
                '   -p     {sample} && \n\n'
                'sort -k 1,1 {go} | join {projpath}/14.enrich/{abbr}.miRNA-geneid - | sed \'s/ /\\t/g\''
                '>{projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt && \\\n'
                'mkdir -p {projpath}/14.enrich/GO_classify && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/go_classification_v3.pl \\\n'
                '   {projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt \\\n'
                '   {projpath}/14.enrich/GO_classify && \n\n'
                'mkdir -p {projpath}/14.enrich/Pathway && \\\n'
                'export PATH={R_v2153}:$PATH && \\\n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \\\n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/runKEGG.pl \\\n'
                '   -directory-gtf {gtf} \\\n'
                '   -gene-id {projpath}/14.enrich/{abbr}.miRNA-geneid \\\n '
                '   -exist-genome y \\\n'
                '   -directory-kobas-output {projpath}/14.enrich/Pathway \\\n'
                '   -directory-genome {refer} \\\n'
                '   -species {keggspe} \\\n'
                '   -samplename {sample} && \n\n'
                'cd {projpath}/14.enrich/Pathway && \\\n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table {projpath}/14.enrich/Pathway/add.${sample}.identify.xls \\\n'
                '   --abbr {keggspe} && \n\n'
                'date +"%D %T ->Finish "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,gtf = self.gtf,sample = sample,abbr = self.prefix,refer = self.refer,
                 pythonExec = self.python2710, R_v2153 = self.R_v2153,LIBRARY_PATH = self.LIBRARY_PATH,
                 PyPakage = self.PyPakage, kobas_v1 = self.kobas_v1, kobas_v2 = self.kobas_v2, kobas =self.kobas,
                 keggspe = self.keggspe)
        
        shell_path = os.path.join(self.projpath,'14.enrich','{}_GO_and_kegg_enrich_animal_single.sh'.format(sample))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_refaimal','diffAnalysis'])
        self.job.add_order(job_name, after_jobs='srna_result')
        

    def GO_and_kegg_enrich_norefanimal(self,each_compare):
        
        code = (
                'date +"%D %T ->Start GO and kegg enrichment " && \\\n\n'
                'cd {projpath}/14.enrich/{each_compare} && \\\n'
                'sort -k 1,1 {projpath}/13.diff/diffAnalysisResult/{each_compare}/{each_compare}.DElist.txt \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID && \n\n'
                'sort -k 1,1 {projpath}/12.target/{abbr}_targets_gene.pairs | \\\n'
                '   join {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNAID - | sort -u \\\n'
                '   > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs && \n\n'
                'awk \'{{print $2}}\' {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-gene.pairs |sort -u \\\n'
                '   >{projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph_v3.pl \\\n'
                '   -i {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -goann {go} \\\n'
                '   -length {length} \\\n'
                '   -o {projpath}/14.enrich/{each_compare}/GO2 \\\n'
                '   -p {each_compare} && \n\n'
                'mkdir -p {projpath}/14.enrich/GO_classify && \\\n'
                'sort -k 1,1 {go} | \\\n'
                '   join {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid - | \\\n'
                '   sed \'s/ /\\t/g\' > {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid.go.txt && \n'
                '{perlExec} {BioScripts}/11Enrichment/go_classification_v3.pl \\\n'
                '   {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid.go.txt \\\n'
                '   {projpath}/14.enrich/each_compare/GO_classify && \n\n'
                'mkdir -p {projpath}/14.enrich/{each_compare}/Pathway && \n'
                'export PATH={R_v2153}:$PATH && \n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/runKEGG_enrich_v1.1.pl  \\\n'
                '   -diff {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -ko   {koann} \\\n'
                '   -g    {each_compare} \\\n'
                '   -t    animal  && \n\n'
                'cd {projpath}/14.enrich/{each_compare}/Pathway && \n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table {projpath}/14.enrich/{each_compare}/{each_compare}.DEG_KEGG_pathway_enrichment_add.xls \\\n' 
                '   --abbr  ko && \\\n\n'
                'date +"%D %T ->Finish GO and kegg enrichment "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,each_compare = each_compare,length = self.length,koann = self.koann,
                 abbr = self.prefix,R_v2153 = self.R_v2153,kobas_v1 = self.kobas_v1,
                 kobas_v2 = self.kobas_v2, LIBRARY_PATH = self.LIBRARY_PATH, kobas = self.kobas,
                 PyPakage = self.PyPakage, pythonExec = self.python2710)

        shell_path = os.path.join(self.projpath,'14.enrich',each_compare,'{}_GO_and_kegg_enrich.sh'.format(each_compare))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_norefanimal','diffAnalysis'])
        self.job.add_order(job_name, after_jobs='srna_result')
        

    def GO_and_kegg_enrich_norefanimal_single(self,sample):
        
        code = (
                'date +"%D %T ->Start single sample GO and kegg enrich analysis " && \\\n\n'
                'cd {projpath}/14.enrich && \\\n'
                'awk \'{{print $2}}\' {projpath}/12.target/{abbr}_targets_gene.pairs | sort -u \\\n' 
                '     >{projpath}/14.enrich/{abbr}.miRNA-geneid && \n\n'
                '{perlExec}  {BioScripts}/11Enrichment/goseq_graph.pl \\\n'
                '   -i      {projpath}/14.enrich/{abbr}.miRNA-geneid \\\n'
                '   -goann  {go} \\\n'
                '   -length {length} \\\n'
                '   -o      {projpath}/14.enrich/GO_enrichment \\\n'
                '   -p      {sample} && \n\n'
                'sort -k 1,1 {go} | join {projpath}/14.enrich/{abbr}.miRNA-geneid - | sed \'s/ /\\t/g\''
                '>{projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt && \\\n'
                'mkdir -p {projpath}/14.enrich/GO_classify && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/go_classification_v3.pl \\\n'
                '   {projpath}/14.enrich/{abbr}.miRNA-geneid.go.txt \\\n'
                '   {projpath}/14.enrich/GO_classify && \n\n'
                'mkdir -p {projpath}/14.enrich/Pathway && \\\n'
                'export PATH={R_v2153}:$PATH && \\\n'
                'export LD_LIBRARY_PATH={LIBRARY_PATH}/R_v2153:$LD_LIBRARY_PATH && \\\n'
                'export PYTHONPATH={PyPakage}:{kobas_v2}/src:{kobas_v1}/src:{kobas}:$PYTHONPATH && \\\n'
                '{perlExec} {BioScripts}/11Enrichment/runKEGG_enrich_v1.1.pl \\\n'               
                '   -diff {projpath}/14.enrich/{each_compare}/{each_compare}.diffmiRNA-geneid \\\n'
                '   -ko   {koann} \\\n'
                '   -g    {sample} \\\n'
                '   -t     animal  && \n\n'
                'cd {projpath}/14.enrich/Pathway && \\\n'
                '{pythonExec} {BioScripts}/11Enrichment/pathway_annotation_flow_parallel_simple_tolerant.py \\\n'
                '   --table {projpath}/14.enrich/Pathway/add.{sample}.identify.xls \\\n'
                '   --abbr ko && \n\n'
                'date +"%D %T ->Finish "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 go = self.go,length = self.length,sample = sample,abbr = self.prefix,koann = self.koann,
                 pythonExec = self.python2710, R_v2153 = self.R_v2153,LIBRARY_PATH = self.LIBRARY_PATH,
                 PyPakage = self.PyPakage, kobas_v1 = self.kobas_v1, kobas_v2 = self.kobas_v2, kobas =self.kobas)

        shell_path = os.path.join(self.projpath,'14.enrich','{}_GO_and_kegg_enrich_norefanimal_single.sh'.format(sample))
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'GO_and_kegg_enrich',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['mirna_target_norefanimal','diffAnalysis'])
        self.job.add_order(job_name, after_jobs='srna_result')


    def srna_result(self):
        
        code = (
                'date +"%D %T ->Start generate srna result " && \\\n'
                '{pythonExec} {Sshow}/srna_result.py \\\n'
                '   --projpath {projpath} \\\n'
                '   --project {project} \\\n'
                '   --sample {sample} \\\n'
                '   --common {common} \\\n'
                '   --org {org} \\\n'
                '   --group {group} \\\n'
                '   --groupname  {groupname} \\\n'
                '   --compare {compare} \\\n'
                '   --go {go} \\\n'
                '   --geneAnn {geneAnn} \\\n'
        ).format(pythonExec = self.python276, Sshow = self.Sshow, project = self.project,
                 projpath = self.projpath, org =self.org, sample = self.sample,common = self.common,
                 groupname = self.groupname, compare = self.compare, group = self.group,
                 go = self.go, geneAnn = self.geneAnn)
        
        if 'noref' in self.org:
            pass
        else:
            if 'refanimal' in self.org:
                code += ('   --utr3 {utr3} '.format(utr3 = self.utr3))
            elif 'refplant' in self.org:
                code += ('   --gene {gene} '.format(gene = self.gene))

        if self.venn_cluster :
            code +=  ('\\\n   --venn_cluster {venn_cluster} && \\\n'
                      'date +"%D %T ->Finish generate srna result "\n').format(venn_cluster = self.venn_cluster)
        else:
            code += ('&& \\\ndate +"%D %T ->Finish generate srna result "\n')
         
        shell_path = os.path.join(self.projpath,self.project+"_sRNA_result",'srna_result.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'srna_result',sched=self.sched)
        self.job.add_atom_job(atom_job)
       

    def srna_report(self):
        
        contract = ''.join(map(self.bracketEscape,list(self.contract)))
        code = (
                'date +"%D %T ->Start generate srna report " && \\\n'
                '{pythonExec} {Sshow}/srna_report.py \\\n'
                '   --projpath {projpath} \\\n'
                '   --project  {project} \\\n'
                '   --contract {contract} \\\n'
                '   --sample {sample} \\\n'
                '   --group {group} \\\n'
                '   --groupname  {groupname} \\\n'
                '   --compare {compare} \\\n'
                '   --common {common} \\\n'
                '   --org {org} '
               ).format(pythonExec = self.python276, Sshow = self.Sshow, project = self.project,
                        projpath = self.projpath,org =self.org, sample = self.sample,
                        common = self.common,contract = contract,groupname = self.groupname,
                        compare = self.compare, group = self.group)

        if self.venn_cluster :
            code +=  ('\\\n   --venn_cluster {venn_cluster} && \\\n'
                      'date +"%D %T ->Finish generate srna report "\n').format(venn_cluster = self.venn_cluster)
        else:
            code += ('&& \\\ndate +"%D %T ->Finish generate srna report "\n')

        shell_path = os.path.join(self.projpath,self.project+"_sRNA_result",'srna_report.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'srna_report',sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='srna_result')


    @staticmethod
    def bracketEscape(character):
        if character in ('(',')','（','）'):
            return '\\'+character
        else:
            return character
