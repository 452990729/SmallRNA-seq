#!/usr/bin/env python
# _*_ coding:utf-8 _*_
 

import os
import glob
from Settings import BioConfig
from Stools import CheckTools, AtomJob, Job, write_shell ,make_dir
from collections import defaultdict


class sRNAClassify(object):

    def __init__(self,args,job):

        self.args         = args
        self.job          = job
        self.projpath     = args.get('projpath') or os.getcwd()
        self.org          = args.get('org').strip()
        self.sched        = args.get('sched')
        self.abbr         = args.get('abbr')
        self.mdspe        = args.get('mdspe') if args.get('mdspe') else ''
        self.NAT          = args.get('NAT')
        self.mode         = args.get('mode')
        self.prefix       = args.get('abbr').strip().split(',')[0]  #NAT 预测的前缀名
        self.sample_list  = CheckTools(self.args).samples
        self.BioScripts   = BioConfig.config.get('srnaenv','BioScripts')
        self.perlExec     = BioConfig.config.get('srnaenv','perl_v5240')
        self.bioperl     = BioConfig.config.get('srnaenv','perl_bio')
        self.miRBase21    = BioConfig.config.get('database','miRBase21')
        self.mirdeep2     = BioConfig.config.get('software','mirdeep2')
        self.srnatoolscli = BioConfig.config.get('software','srnatoolscli')
        self.ViennaRNA2   = BioConfig.config.get('software','ViennaRNA2')
        self.bowtie1      = BioConfig.config.get('software','bowtie1')
#        self.miREvo       = BioConfig.config.get('software','miREvo')
#        self.known_TAS_db = BioConfig.config.get('database','known_TAS')   
#        self.rRNA         = CheckTools(self.args).getRefConfig().get('rRNA')
#        self.tRNA         = CheckTools(self.args).getRefConfig().get('tRNA')
#        self.snRNA        = CheckTools(self.args).getRefConfig().get('snRNA')
#        self.snoRNA       = CheckTools(self.args).getRefConfig().get('snoRNA')
#        self.repdir       = CheckTools(self.args).getRefConfig().get('repdir')
#        self.refer        = CheckTools(self.args).getRefConfig().get('refer')        
#        self.gtf          = CheckTools(self.args).getRefConfig().get('gtf')
#        self.exon         = CheckTools(self.args).getRefConfig().get('exon')
#        self.intron       = CheckTools(self.args).getRefConfig().get('intron')


    def run(self):

        self.get_hairpin_mature_seq()
        self.known_miRNA_quantify()
        self.known_miRNA_csv2pdf()            
        self.known_miRNA_stat()
        '''
        self.ncRNA_map_and_stat()
        if 'noref' in self.org.lower():
            self.predict_novel_miRNA()           
        else:
            self.repeat_mapping()
            if 'plant' in self.org.lower():
                self.NAT_siRNA_plant()
            if 'hsa' not in self.mdspe.lower() and 'hsa' not in self.abbr:
                self.exon_and_intron_build_index()
            self.gene_map_and_stat()
            self.predict_novel_miRNA()
        self.novel_miRNA_quantify()
        self.novel_miRNA_stat()
        if self.org.lower() == 'refplant':
            self.known_TAS_and_novel_TAS_analysis()
            self.TAS_mapping_and_TAS_stat()
        self.sRNA_Category()
        '''

    def get_hairpin_mature_seq(self):

        miRNA_Reference_species_path  = os.path.join(self.projpath,'3.known','miRbase')
        make_dir(miRNA_Reference_species_path)
        miRNA_Reference_species_list  = self.abbr.strip().split(',')
        with open(os.path.join(miRNA_Reference_species_path,'mir_list.txt'),'w+') as file:
            file.writelines([line +'\n' for line in miRNA_Reference_species_list]) 

        code = (
                'date +"%D %T -> Start get_hairpin_mature_seq " && \\\n'
                '{perlExec} {BioScripts}/3known/get_hairpin_or_mature_from_miRBase.pl \\\n'
                '    {projpath}/3.known/miRbase/mir_list.txt \\\n'
                '    {miRBase21}/hairpin.fa \\\n'
                '    > {projpath}/3.known/miRbase/ref_hairpin.fa && \n'
                '{perlExec} {BioScripts}/3known/get_hairpin_or_mature_from_miRBase.pl \\\n'
                '    {projpath}/3.known/miRbase/mir_list.txt \\\n'
                '    {miRBase21}/mature.fa \\\n'
                '     > {projpath}/3.known/miRbase/ref_mature.fa && \n\n'
                '{perlExec} {BioScripts}/3known/known_hairpin_mature_rmdup.pl \\\n'
                '    {projpath}/3.known/miRbase/ref_hairpin.fa \\\n'
                '    {projpath}/3.known/miRbase/ref_mature.fa \\\n'
                '    {projpath}/3.known/miRbase/known_hairpin_ref.fa \\\n'
                '    {projpath}/3.known/miRbase/known_mature_ref.fa && \n\n'
                'cat {projpath}/2.map/*/*.map.collapse.fa > {projpath}/3.known/miRbase/all.sample.map.collapse.fa && \n'
                'date +"%D %T -> Finish get_hairpin_mature_seq "\n'
            ).format(projpath  = self.projpath,BioScripts = self.BioScripts,perlExec = self.bioperl,
                     miRBase21 = self.miRBase21)
        shell_path = os.path.join(self.projpath,'3.known','miRbase','get_hairpin_mature_seq.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'get_hairpin_mature_seq', sched=self.sched)
        self.job.add_atom_job(atom_job)
        #self.job.add_order(job_name, previous_jobs='qc_report')


    def known_miRNA_quantify(self):
        code = (
                'date +"%D %T -> Start known miRNA quantify " && \\\n'
                'cd {projpath}/3.known && \n'
                'export PATH={mirdeep2}:$PATH \n'
                '{perlExec} {BioScripts}/3known/quantifier_gb_v2.pl \\\n'
                '   -p {projpath}/3.known/miRbase/known_hairpin_ref.fa \\\n'
                '   -m {projpath}/3.known/miRbase/known_mature_ref.fa \\\n'
                '   -r {projpath}/3.known/miRbase/all.sample.map.collapse.fa \\\n'
                '   -y known_miRNAs_expressed \\\n'
                '   -g 0 \\\n'
                '   -T 5 && \\\n'
                'date +"%D %T -> Finish known miRNA quantify "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 mirdeep2 = self.mirdeep2)
        shell_path = os.path.join(self.projpath,'3.known','known_miRNA_quantify.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'known_miRNA_quantify', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='get_hairpin_mature_seq')


    def known_miRNA_csv2pdf(self):
        code = (
                'date +"%D %T ->Start known miRNA csv2pdf " && \\\n'
                'export PATH={srnatoolscli}:{ViennaRNA2}/bin:$PATH \n'
                '{perlExec} {BioScripts}/3known/csv2pdf_gb.pl \\\n'
                '   {projpath}/3.known/known_miRNAs_expressed.csv \\\n'
                '   {projpath}/3.known/miRbase/known_hairpin_ref.fa \\\n'
                '   {projpath}/3.known/miRbase/known_mature_ref.fa \\\n'
                '   {projpath}/3.known/known_miRNAs && \n\n'
                'awk -F "\\t" -v OFS="\\t" \'{{if($4>0){{print $1,$3}}}}\' {projpath}/3.known/known_miRNAs_expressed.csv '
                '> {projpath}/3.known/known_miRNAs/hairpin_mature.pairs && \n\n'
                'head -1 {projpath}/3.known/known_miRNAs/miRNAs_expressed_known.csv | awk \'{{num=(NF-4)/2;printf("miRNA");'
                'for(i=1;i<=num;i++){{printf("\\t"$(i+4))}};printf("\\n");}}\' > {projpath}/3.known/known_miRNAs/mature.readcount && \n\n'
                'awk \'{{if(NR>1){{num=(NF-4)/2;printf($1"\\t"$2);for(i=1;i<=num;i++){{printf("\\t"$(i+4))}};printf("\\n")}}}}\' '
                '{projpath}/3.known/known_miRNAs/miRNAs_expressed_known.csv | sort -k 1,1 -k 2nr,2 | '
                'awk \'{{if($1!=name){{num=NF-2;printf($1);for(i=1;i<=num;i++){{printf("\\t"$(i+2));}}printf("\\n");name=$1}}}}\' '
                '>>{projpath}/3.known/known_miRNAs/mature.readcount  && \n\n'
                'awk \'{{if(/^>/){{title=$1}}else{{if(title){{if(/total read count/ && $(NF)>0)'
                '{{print title;title="";marker=1}}else{{marker=0}}}}if(marker){{print}}}}}}\' '
                '{projpath}/3.known/expression_analyses/miRBase.mrd >{projpath}/3.known/known_miRNAs/miRBase.mrd && \n\n'
                'date +"%D %T ->Finish known miRNA csv2pdf"\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 srnatoolscli = self.srnatoolscli,ViennaRNA2 = self.ViennaRNA2)
        shell_path = os.path.join(self.projpath,'3.known','known_miRNA_csv2pdf.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]   
        atom_job = AtomJob(job_name, shell_path,1,'known_miRNA_csv2pdf', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='known_miRNA_quantify')


    def known_miRNA_stat(self):
        code = (
                'date +"%D %T ->Start known miRNA stat " && \\\n\n'
                'cd {projpath}/3.known && \\\n'
                '{perlExec} {BioScripts}/3known/genebwt12count.pl \\\n'
                '   -i {projpath}/3.known/expression_analyses/*_mapped.bwt.ka \\\n'
                '   -r {projpath}/3.known/expression_analyses/precursor.converted \\\n'
                '   -t  known_miRNA \\\n'
                '   -o  {projpath}/3.known/known_miRNAs \\\n'
                '   -u -s -W && \n\n'
                'awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped mature\\t"$0}}}} \' \\\n'
                '   {projpath}/3.known/known_miRNAs/known_miRNA.mapmat.stat \\\n'
                '   > {projpath}/3.known/known_miRNAs/known_miRNA.map.stat && \n\n'
                '{perlExec} {BioScripts}/3known/stat_known_miRNA_pre.pl \\\n'
                '   {projpath}/3.known/known_miRNAs/miRBase.mrd \\\n'
                '   >{projpath}/3.known/known_miRNAs/known_miRNA.mapref.stat  && \n'
                'awk \'{{if(NR==2){{print "Mapped hairpin\\t"$0}}}}\' \\\n'
                '   {projpath}/3.known/known_miRNAs/known_miRNA.mapref.stat \\\n' 
                '   >> {projpath}/3.known/known_miRNAs/known_miRNA.map.stat && \n\n'
                'awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=int($i+0.5)}}printf("Mapped uniq sRNA\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"int($i+0.5))}}printf("\\n");}}}}\' \\\n'
                '   {projpath}/3.known/known_miRNAs/known_miRNA.uc.stat \\\n'
                '   >>{projpath}/3.known/known_miRNAs/known_miRNA.map.stat && \n\n'
                'awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=int($i+0.5)}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"int($i+0.5))}}printf("\\n");}}}}\' \\\n' 
                '   {projpath}/3.known/known_miRNAs/known_miRNA.rc.stat \\\n'
                '   >>{projpath}/3.known/known_miRNAs/known_miRNA.map.stat && \n\n'
                'ln -sf {projpath}/3.known/expression_analyses/*.unmap.fas \\\n'
                '   {projpath}/3.known/known_miRNAs/known_miRNA.unmap.fas && \n\n'
                '{perlExec} {BioScripts}/3known/bwt12collapse.pl \\\n'
                '   {projpath}/3.known/expression_analyses/*_mapped.bwt.k1 \\\n'
                '   >{projpath}/3.known/known_miRNAs//known_miRNA.map.fas && \n\n'
                '{perlExec} {BioScripts}/3known/bwt_base_bias_plot.pl \\\n'
                '    -i {projpath}/3.known/expression_analyses/*_mapped.bwt.k1 && \n\n'
                'date +"%D %T ->Finish known miRNA stat "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec)
        shell_path = os.path.join(self.projpath,'3.known','known_miRNA_stat.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'known_miRNA_stat', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='known_miRNA_csv2pdf')


    def ncRNA_map_and_stat(self):
        input_dir  = os.path.join(self.projpath,'4.ncRNA/input')
        output_dir = os.path.join(self.projpath,'4.ncRNA/output')
        make_dir(input_dir)
        make_dir(output_dir)
        code = (
                'date +"%D %T ->Stat ncRNA map and stat " && \\\n'
                'echo ================== Run rRNA mapping and Stat ==================\n'
                'ln -sf {rRNA} \\\n'
                '       {input_dir}/rRNA.fna && \n'
                '{bowtie1}/bowtie-build \\\n'
                '       {input_dir}/rRNA.fna \\\n'
                '       {input_dir}/rRNA.fna && \n'                
                '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                '   -f {input_dir}/rRNA.fna \\\n'
                '      {projpath}/3.known/known_miRNAs/known_miRNA.unmap.fas \\\n'
                '      {output_dir}/rRNA.bwt \\\n'
                '   --un {output_dir}/rRNA.unmap.fas \\\n'
                '   2>{output_dir}/rRNA.mapping.log && \n\n'
                'if [[ ! -z {output_dir}/rRNA.bwt  ]]; then \n'
                '   {perlExec} {BioScripts}/4ncRNA/bwt12collapse.pl \\\n'
                '       {output_dir}/rRNA.bwt \\\n'
                '       >{output_dir}/rRNA.map.fas && \n'
                '   cat {output_dir}/rRNA.map.fas \\\n'
                '       > {output_dir}/ncRNA.map.fas && \n'
                '   {perlExec}  {BioScripts}/4ncRNA/genebwt12count.pl \\\n'
                '       -i {output_dir}/rRNA.bwt \\\n'
                '       -r {input_dir}/rRNA.fna \\\n'
                '       -t rRNA \\\n'
                '       -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/rRNA.mapref.stat \\\n'
                '       >{output_dir}/rRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/rRNA.uc.stat \\\n'
                '       >>{output_dir}/rRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/rRNA.rc.stat \\\n'
                '       >>{output_dir}/rRNA.map.stat \n'
                'fi && \n\n'
                'echo ================== Run tRNA mapping and Stat ==================\n'
                'ln -sf {tRNA} \\\n'
                '       {input_dir}/tRNA.fna && \n'
                '{bowtie1}/bowtie-build \\\n'
                '       {input_dir}/tRNA.fna \\\n'
                '       {input_dir}/tRNA.fna && \n'                
                '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                '   -f {input_dir}/tRNA.fna  \\\n'
                '      {output_dir}/rRNA.unmap.fas \\\n'
                '      {output_dir}/tRNA.bwt \\\n'
                '   --un {output_dir}/tRNA.unmap.fas \\\n'
                '   2>{output_dir}/tRNA.mapping.log && \n\n'
                'if [[ ! -z {output_dir}/tRNA.bwt  ]]; then \n'
                '   {perlExec} {BioScripts}/4ncRNA/bwt12collapse.pl \\\n'
                '       {output_dir}/tRNA.bwt \\\n'
                '       >{output_dir}/tRNA.map.fas && \n'
                '   cat {output_dir}/tRNA.map.fas \\\n'
                '       >> {output_dir}/ncRNA.map.fas && \n'
                '   {perlExec}  {BioScripts}/4ncRNA/genebwt12count.pl \\\n'
                '       -i {output_dir}/tRNA.bwt \\\n'
                '       -r {input_dir}/tRNA.fna \\\n'
                '       -t tRNA \\\n'
                '       -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/tRNA.mapref.stat \\\n'
                '       >{output_dir}/tRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/tRNA.uc.stat \\\n'
                '       >>{output_dir}/tRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/tRNA.rc.stat \\\n'
                '       >>{output_dir}/tRNA.map.stat \n'
                'fi && \n\n'
                'echo ================== Run snRNA mapping and Stat ==================\n'
                'ln -sf {snRNA} \\\n'
                '       {input_dir}/snRNA.fna && \n'
                '{bowtie1}/bowtie-build \\\n'
                '       {input_dir}/snRNA.fna \\\n'
                '       {input_dir}/snRNA.fna && \n'                
                '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                '   -f {input_dir}/snRNA.fna \\\n'
                '      {output_dir}/tRNA.unmap.fas \\\n'
                '      {output_dir}/snRNA.bwt \\\n'
                '   --un {output_dir}/snRNA.unmap.fas \\\n'
                '   2>{output_dir}/snRNA.mapping.log && \n\n'
                'if [[ ! -z {output_dir}/snRNA.bwt  ]]; then \n'
                '   {perlExec} {BioScripts}/4ncRNA/bwt12collapse.pl \\\n'
                '       {output_dir}/snRNA.bwt \\\n'
                '       >{output_dir}/snRNA.map.fas && \n'
                '   cat {output_dir}/snRNA.map.fas \\\n'
                '       >> {output_dir}/ncRNA.map.fas && \n'
                '   {perlExec}  {BioScripts}/4ncRNA/genebwt12count.pl \\\n'
                '       -i {output_dir}/snRNA.bwt \\\n'
                '       -r {input_dir}/snRNA.fna \\\n'
                '       -t snRNA \\\n'
                '       -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/snRNA.mapref.stat \\\n'
                '       >{output_dir}/snRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/snRNA.uc.stat \\\n'
                '       >>{output_dir}/snRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/snRNA.rc.stat \\\n'
                '       >>{output_dir}/snRNA.map.stat \n'
                'fi && \n\n'
                'echo ================== Run snoRNA mapping and Stat ==================\n'
                'ln -sf {snoRNA} \\\n'
                '       {input_dir}/snoRNA.fna && \n'
                '{bowtie1}/bowtie-build \\\n'
                '       {input_dir}/snoRNA.fna \\\n'
                '       {input_dir}/snoRNA.fna && \n'                
                '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                '   -f {input_dir}/snoRNA.fna \\\n'
                '      {output_dir}/snRNA.unmap.fas \\\n'
                '      {output_dir}/snoRNA.bwt \\\n'
                '   --un {output_dir}/snoRNA.unmap.fas \\\n'
                '   2>{output_dir}/snoRNA.mapping.log && \n\n'
                'if [[ ! -z {output_dir}/snoRNA.bwt  ]]; then \n'
                '   {perlExec} {BioScripts}/4ncRNA/bwt12collapse.pl \\\n'
                '       {output_dir}/snoRNA.bwt \\\n'
                '       >{output_dir}/snoRNA.map.fas && \n'
                '   cat {output_dir}/snoRNA.map.fas \\\n'
                '       >> {output_dir}/ncRNA.map.fas && \n'
                '   {perlExec}  {BioScripts}/4ncRNA/genebwt12count.pl \\\n'
                '       -i {output_dir}/snoRNA.bwt \\\n'
                '       -r {input_dir}/snoRNA.fna \\\n'
                '       -t snoRNA \\\n'
                '       -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/snoRNA.mapref.stat \\\n'
                '       >{output_dir}/snoRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/snoRNA.uc.stat \\\n'
                '       >>{output_dir}/snoRNA.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/snoRNA.rc.stat \\\n'
                '       >>{output_dir}/snoRNA.map.stat  \n'
                'fi && \n\n'
                '{perlExec} {BioScripts}/4ncRNA/paste_col.pl \\\n'
                '       {output_dir}/rRNA.uc.stat  \\\n'
                '       {output_dir}/tRNA.uc.stat \\\n'
                '       {output_dir}/snRNA.uc.stat \\\n'
                '       {output_dir}/snoRNA.uc.stat \\\n'
                '        >{output_dir}/uc.stat && \n'
                '{perlExec} {BioScripts}/4ncRNA/paste_col.pl \\\n'
                '       {output_dir}/rRNA.rc.stat  \\\n'
                '       {output_dir}/tRNA.rc.stat \\\n'
                '       {output_dir}/snRNA.rc.stat \\\n'
                '       {output_dir}/snoRNA.rc.stat \\\n'
                '        >{output_dir}/rc.stat && \n'
                'date +"%D %T ->Finish ncRNA map and stat "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 input_dir = input_dir,output_dir = output_dir,rRNA = self.rRNA,tRNA = self.tRNA,
                 snRNA = self.snRNA, snoRNA = self.snoRNA, bowtie1 = self.bowtie1)
        shell_path = os.path.join(self.projpath,'4.ncRNA','ncRNA_map_and_stat.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'ncRNA_map_and_stat', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='known_miRNA_stat')

    def repeat_mapping(self):
        code = (
                'date +"%D %T ->Stat repeat_mapping " && \\\n'
                '{perlExec} {BioScripts}/5repeat/repeat_map.pl \\\n'
                '   {projpath}/4.ncRNA/output/snoRNA.unmap.fas \\\n'
                '   {repdir} \\\n'
                '   {projpath}/5.repeat/repeat_result && \n'
                'date +"%D %T ->Finish repeat_mapping "\n'
    ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
             repdir = self.repdir)
        shell_path = os.path.join(self.projpath,'5.repeat','repeat_mapping.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'repeat_mapping', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='ncRNA_map_and_stat')

    def NAT_siRNA_plant(self):
       
        shell_path = os.path.join(self.projpath,'6.NAT','NAT_siRNA_plant.sh')
        if self.NAT.lower() != 'other':       
            code  = (
                    'date +"%D %T ->Start Natural Antisense Transcripts Screening" && \\\n'
                    '{perlExec} {BioScripts}/6plant_NAT/NAT_map.pl \\\n'
                    '   -i {projpath}/5.repeat/repeat_result/repeat.unmap.fas \\\n'
                    '   -s {NAT} \\\n'
                    '   -o {projpath}/6.NAT && \n'
                    'date +"%D %T ->Finish Natural Antisense Transcripts Screening "\n'
            ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                     NAT = self.NAT)
            write_shell(code,shell_path)
        else:
            '''
                NAT_predict_map.pl生成运行的shell脚本
            '''
            cmd  = (
                    '{perlExec} {BioScripts}/6plant_NAT/NAT_predict_map.pl \\\n'
                    '   -q   {projpath}/5.repeat/repeat_result/repeat.unmap.fas \\\n'
                    '   -o   {prefix} \\\n'
                    '   -gtf {gtf} \\\n'
                    '   -genome {refer} \\\n'
                    '   -outdir {projpath}/6.NAT \\\n'
                    '   > {projpath}/6.NAT/NAT_siRNA_plant.sh'
            ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                     prefix = self.prefix,gtf = self.gtf,refer = self.refer)
            assert not os.system(cmd)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'NAT_siRNA_plant', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='repeat_mapping')
    
    def exon_and_intron_build_index(self):

        input_dir  = os.path.join(self.projpath,'7.gene/input')        
        output_dir = os.path.join(self.projpath,'7.gene/output')
        make_dir(input_dir)
        make_dir(output_dir)

        code = (
                'date +"%D %T ->Start exon and intron build index" && \\\n'
                'ln -sf {exon} \\\n'
                '       {input_dir}/exon.fna && \n'
                '{bowtie1}/bowtie-build \\\n'
                '       {input_dir}/exon.fna \\\n'
                '       {input_dir}/exon.fna && \n'                
                'ln -sf {intron} \\\n'
                '       {input_dir}/intron.fna && \n'
                '{bowtie1}/bowtie-build \\\n'
                '       {input_dir}/intron.fna \\\n'
                '       {input_dir}/intron.fna && \n'                
                'date +"%D %T ->Finish exon and intron build index "\n'
        ).format(input_dir = input_dir, output_dir = output_dir, bowtie1 = self.bowtie1,
                 exon = self.exon, intron = self.intron)
        shell_path = os.path.join(self.projpath,'7.gene','exon_and_intron_build_index.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name,shell_path,1,'exon_and_intron_build_index', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, after_jobs='gene_map_and_stat')        

            
    def gene_map_and_stat(self):

        input_dir  = os.path.join(self.projpath,'7.gene/input')
        output_dir = os.path.join(self.projpath,'7.gene/output')
        make_dir(input_dir)
        make_dir(output_dir)

        if 'noref' in self.org.lower():
            unmapping_fas = '{projpath}/4.ncRNA/output/snoRNA.unmap.fas'.format(projpath = self.projpath)
            previous_jobs = 'ncRNA_map_and_stat'
        else:            
            if 'plant' in self.org.lower():
                unmapping_fas = '{projpath}/6.NAT/output/NAT.unmap.fas'.format(projpath = self.projpath)
                previous_jobs = 'NAT_siRNA_plant'
            else:
                unmapping_fas = '{projpath}/5.repeat/repeat_result/repeat.unmap.fas'.format(projpath = self.projpath)
                previous_jobs = 'repeat_mapping'

            if 'hsa' not in self.mdspe.lower() and 'hsa' not in self.abbr :
                exon_fna      = '{input_dir}/exon.fna'.format(input_dir = input_dir) 
                intron_fna    = '{input_dir}/intron.fna'.format(input_dir = input_dir)
            else:
                exon_fna      = self.exon
                intron_fna    = self.intron

        code = (
                'date +"%D %T -> Start gene mapping and stat" && \\\n'
                'echo ================== Run exon mapping and stat ==================\n'
                '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                '   -f {exon_fna} \\\n'
                '      {unmapping_fas} \\\n'
                '      {output_dir}/exon.bwt \\\n'
                '   --un {output_dir}/exon.unmap.fas \\\n'
                '   2>{output_dir}/exon.mapping.log && \n\n'
                'if [[ ! -z {output_dir}/exon.bwt  ]]; then \n'
                '   {perlExec} {BioScripts}/4ncRNA/bwt12collapse.pl \\\n'
                '       {output_dir}/exon.bwt \\\n'
                '       >{output_dir}/exon.map.fas && \n'
                '   cat {output_dir}/exon.map.fas \\\n'
                '       > {output_dir}/gene.map.fas && \n'
                '   {perlExec}  {BioScripts}/4ncRNA/genebwt12count.pl \\\n'
                '       -i {output_dir}/exon.bwt \\\n'
                '       -r {exon_fna} \\\n'
                '       -t exon \\\n'
                '       -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/exon.mapref.stat \\\n'
                '       >{output_dir}/exon.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/exon.uc.stat \\\n'
                '       >>{output_dir}/exon.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/exon.rc.stat \\\n'
                '       >>{output_dir}/exon.map.stat \n'
                'fi && \n\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 output_dir = output_dir,unmapping_fas = unmapping_fas,exon_fna = exon_fna,
                 bowtie1 = self.bowtie1)
                 
        if 'hsa' in self.abbr:
            
            intron_prefix = os.path.splitext(intron_fna)[0]
            
            code += (
                    'echo ================== Run intron mapping and Stat ==================\n'
                    '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                    '   -f {intron_prefix}_1.fa  \\\n'
                    '      {output_dir}/exon.unmap.fas \\\n'
                    '      {output_dir}/intron1.bwt \\\n'
                    '   --un {output_dir}/intron1.unmap.fas \\\n'
                    '   2>{output_dir}/intron1.mapping.log && \n'
                    '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                    '   -f {intron_prefix}_2.fa  \\\n'
                    '      {output_dir}/intron1.unmap.fas \\\n'
                    '      {output_dir}/intron2.bwt \\\n'
                    '   --un {output_dir}/intron2.unmap.fas \\\n'
                    '   2>{output_dir}/intron2.mapping.log && \n'                
                    '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                    '   -f {intron_prefix}_3.fa  \\\n'
                    '      {output_dir}/intron2.unmap.fas \\\n'
                    '      {output_dir}/intron3.bwt \\\n'
                    '   --un {output_dir}/intron3.unmap.fas \\\n'
                    '   2>{output_dir}/intron3.mapping.log && \n'
                    '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                    '   -f {intron_prefix}_4.fa  \\\n'
                    '      {output_dir}/intron3.unmap.fas \\\n'
                    '      {output_dir}/intron4.bwt \\\n'
                    '   --un {output_dir}/intron4.unmap.fas \\\n'
                    '   2>{output_dir}/intron4.mapping.log && \n'
                    '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                    '   -f {intron_prefix}_5.fa  \\\n'
                    '      {output_dir}/intron4.unmap.fas \\\n'
                    '      {output_dir}/intron5.bwt \\\n'
                    '   --un {output_dir}/intron5.unmap.fas \\\n'
                    '   2>{output_dir}/intron5.mapping.log && \n'
                    'cat {projpath}/7.gene/output/intron*.bwt \\\n'
                    '   >{projpath}/7.gene/output/intron.bwt && \n'
                    'ln -sf {output_dir}/intron5.unmap.fas {output_dir}/intron.unmap.fas && \n\n'
            ).format(bowtie1 = self.bowtie1,intron_fna = intron_fna,output_dir = output_dir,
                     projpath = self.projpath,intron_prefix = intron_prefix)

        else:
            code += (
                    'echo ================== Run intron mapping and Stat ==================\n'
                    '{bowtie1}/bowtie -p 10 -v 0 -k 1 \\\n'
                    '   -f {intron_fna}  \\\n'
                    '      {output_dir}/exon.unmap.fas \\\n'
                    '      {output_dir}/intron.bwt \\\n'
                    '   --un {output_dir}/intron.unmap.fas \\\n'
                    '   2>{output_dir}/intron.mapping.log && \n\n'
            ).format(bowtie1 = self.bowtie1,intron_fna = intron_fna,output_dir = output_dir)
        
        
        code +=(
                'if [[ ! -z {output_dir}/intron.bwt  ]]; then \n'
                '   {perlExec} {BioScripts}/4ncRNA/bwt12collapse.pl \\\n'
                '       {output_dir}/intron.bwt \\\n'
                '       >{output_dir}/intron.map.fas && \n'
                '   cat {output_dir}/intron.map.fas \\\n'
                '       >> {output_dir}/gene.map.fas && \n'
                '   {perlExec}  {BioScripts}/4ncRNA/genebwt12count.pl \\\n'
                '       -i {output_dir}/intron.bwt \\\n'
                '       -r {intron_fna} \\\n'
                '       -t intron \\\n'
                '       -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/intron.mapref.stat \\\n'
                '       >{output_dir}/intron.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/intron.uc.stat \\\n'
                '       >>{output_dir}/intron.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/intron.rc.stat \\\n'
                '       >>{output_dir}/intron.map.stat \n'
                'fi && \n\n'
                'sed -n \'2,$p\' {output_dir}/intron.uc.stat | cat {output_dir}/exon.uc.stat - \\\n'
                '       > {output_dir}/uc.stat && \n'
                'sed -n \'2,$p\' {output_dir}/intron.rc.stat | cat {output_dir}/exon.rc.stat - \\\n'
                '       > {output_dir}/rc.stat && \n\n'
                'date +"%D %T -> Finish gene mapping and stat " \n'
        ).format(BioScripts = self.BioScripts,perlExec = self.perlExec,output_dir = output_dir,intron_fna = intron_fna)
        
        shell_path = os.path.join(self.projpath,'7.gene','gene_map_and_stat.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'gene_map_and_stat', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=previous_jobs)
        
        
        
    

    def predict_novel_miRNA(self):

        '''
            k pipline_predict.sh parameter 
            a read is allowed to map up to this number of positions in the genome, 
            default=5; for plant, 15 is recommaned
        '''        

        miRNA_Reference_species_path  = os.path.join(self.projpath,'8.novel','miRbase')
        make_dir(miRNA_Reference_species_path)
        genome = os.path.basename(self.refer)
 
        if 'plant' in self.org.lower():
            mature_seq   = '{miRBase21}/mature_plant.fa'.format(miRBase21 = self.miRBase21)
            k = 15
        else:
            mature_seq   = '{miRBase21}/mature_animal.fa'.format(miRBase21 = self.miRBase21)
            k = 5

        if 'noref' in self.org.lower():
            unmap_fas    = '{projpath}/4.ncRNA/output/snoRNA.unmap.fas'.format(projpath = self.projpath)
            previous_job = 'ncRNA_map_and_stat'
        else:   
            unmap_fas    = '{projpath}/7.gene/output/intron.unmap.fas'.format(projpath = self.projpath)
            previous_job = 'gene_map_and_stat'

        code = (
                'date +"%D %T ->Start predict novel miRNA" && \\\n'
                'ln -sf {mature_seq} \\\n'
                '       {miRNA_Reference_species_path}/novel_mature_ref.fa && \n'
                'ln -sf {refer}* \\\n'
                '       {miRNA_Reference_species_path} && \n'
                'export MIREVO={miREvo} \n'
                'export PERL5LIB={miREvo}:$PERL5LIB \n'
                'export PATH={srnatoolscli}:$PATH \n'
                'export PATH={bowtie1}:$PATH \n'
                '{BioScripts}/8novel/pipline_predict.sh \\\n'
                '       -i {unmap_fas} \\\n'
                '       -r {miRNA_Reference_species_path}/{genome} \\\n'
                '       -M {miRNA_Reference_species_path}/novel_mature_ref.fa \\\n'
                '       -m {mode} \\\n'
                '       -k {k} \\\n'
                '       -p 10 \\\n'
                '       -g 50000 \\\n'
                '       -w {abbr} \\\n'
                '       -o {projpath}/8.novel/predict_novel_miRNA && \n'
                'date +"%D %T ->Finish predict novel miRNA "\n'
        ).format(mature_seq = mature_seq,miRNA_Reference_species_path = miRNA_Reference_species_path,refer = self.refer, 
                 miREvo = self.miREvo, srnatoolscli = self.srnatoolscli,bowtie1 = self.bowtie1,unmap_fas = unmap_fas,
                 genome = genome,mode = self.mode,k = k, BioScripts = self.BioScripts,projpath = self.projpath,abbr=self.abbr)
        shell_path = os.path.join(self.projpath,'8.novel','predict_novel_miRNA.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'predict_novel_miRNA', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=previous_job)
                    
    def novel_miRNA_quantify(self):

        if 'noref' in self.org.lower():
            unmap_fas    = '{projpath}/4.ncRNA/output/snoRNA.unmap.fas'.format(projpath = self.projpath)
        else:   
            unmap_fas    = '{projpath}/7.gene/output/intron.unmap.fas'.format(projpath = self.projpath)
        
        code = (
                'date +"%D %T ->Start novel miRNA quantify " && \\\n'
                'cd {projpath}/8.novel && \\\n'
                'export PATH={mirdeep2}:$PATH \n'
                '{perlExec} {BioScripts}/8novel/novel_hp_mat_remove_dup.pl \\\n'
                '   {projpath}/8.novel/predict_novel_miRNA/predict.result.csv \\\n'
                '   {projpath}/8.novel/miRbase/predict_hairpin.fa \\\n' 
                '   {projpath}/8.novel/miRbase/predict_mature.fa \\\n'
                '   {projpath}/8.novel/miRbase/predict_hairpin.pos  && \n'
                '{perlExec} {BioScripts}/8novel/quantifier_gb_v2.pl \\\n'
                '   -p {projpath}/8.novel/miRbase/predict_hairpin.fa \\\n'
                '   -m {projpath}/8.novel/miRbase/predict_mature.fa \\\n'
                '   -r {unmap_fas} \\\n'
                '   -g 0 -T 10 \\\n'
                '   -y novel_miRNAs_expressed && \n'
                'date +"%D %T ->Finish novel miRNA quantify "\n' 
        ).format(mirdeep2 = self.mirdeep2, perlExec = self.perlExec,BioScripts = self.BioScripts,
                 projpath = self.projpath,unmap_fas = unmap_fas) 
        shell_path = os.path.join(self.projpath,'8.novel','novel_miRNA_quantify.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'novel_miRNA_quantify', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='predict_novel_miRNA')


    def novel_miRNA_stat(self):
        
        code = (
                'date +"%D %T ->Start novel miRNA stat " && \\\n\n'
                'cd {projpath}/8.novel && \\\n'
                'export PATH={srnatoolscli}:{ViennaRNA2}/bin:$PATH && \n'
                '{perlExec} {BioScripts}/8novel/cvs2result.pl \\\n'
                '   {projpath}/8.novel/novel_miRNAs_expressed.csv \\\n'
                '   {projpath}/8.novel/miRbase/predict_hairpin.fa \\\n'
                '   {projpath}/8.novel/miRbase/predict_mature.fa \\\n'
                '   {projpath}/8.novel/miRbase/predict_hairpin.pos \\\n'
                '   {projpath}/8.novel/novel_miRNAs && \n\n'
                'awk -F "\\t" -v OFS="\\t" \'{{if($4>0){{print $1,$3}}}}\' {projpath}/8.novel/novel_miRNAs_expressed.csv '
                '> {projpath}/8.novel/novel_miRNAs/hairpin_mature.pairs && \n\n'
                'head -1 {projpath}/8.novel/novel_miRNAs/miRNAs_expressed_novel.csv | awk \'{{num=(NF-4)/2;printf("miRNA");'
                'for(i=1;i<=num;i++){{printf("\\t"$(i+4))}};printf("\\n");}}\' > {projpath}/8.novel/novel_miRNAs/mature.readcount && \n\n'
                'awk \'{{if(NR>1 && $1!~/*$/){{num=(NF-4)/2;printf($1"\\t"$2);for(i=1;i<=num;i++){{printf("\\t"$(i+4))}};printf("\\n")}}}}\' '
                '{projpath}/8.novel/novel_miRNAs/miRNAs_expressed_novel.csv | sort -k 1,1 -k 2nr,2 | '
                'awk \'{{if($1!=name){{num=NF-2;printf($1);for(i=1;i<=num;i++){{printf("\\t"$(i+2));}}printf("\\n");name=$1}}}}\' '
                '>>{projpath}/8.novel/novel_miRNAs/mature.readcount  && \n\n'
                'head -1 {projpath}/8.novel/novel_miRNAs/miRNAs_expressed_novel.csv | awk \'{{num=(NF-4)/2;printf("miRNA");'
                'for(i=1;i<=num;i++){{printf("\\t"$(i+4))}};printf("\\n");}}\' > {projpath}/8.novel/novel_miRNAs/star.readcount && \n\n'
                'awk \'{{if(NR>1 && $1~/*$/){{num=(NF-4)/2;printf($1"\\t"$2);for(i=1;i<=num;i++){{printf("\\t"$(i+4))}};printf("\\n")}}}}\' '
                '{projpath}/8.novel/novel_miRNAs/miRNAs_expressed_novel.csv | sort -k 1,1 -k 2nr,2 | '
                'awk \'{{if($1!=name){{num=NF-2;printf($1);for(i=1;i<=num;i++){{printf("\\t"$(i+2));}}printf("\\n");name=$1}}}}\' '
                '>>{projpath}/8.novel/novel_miRNAs/star.readcount  && \n\n'
                'awk \'{{if(/^>/){{title=$1}}else{{if(title){{if(/total read count/ && $(NF)>0)'
                '{{print title;title="";marker=1}}else{{marker=0}}}}if(marker){{print}}}}}}\' '
                '{projpath}/8.novel/expression_analyses/miRBase.mrd >{projpath}/8.novel/novel_miRNAs/miRBase.mrd && \n\n'
                '{perlExec} {BioScripts}/8novel/genebwt12count.pl \\\n'
                '   -i {projpath}/8.novel/expression_analyses/novel_mapped.bwt.ka \\\n'
                '   -r {projpath}/8.novel/expression_analyses/precursor.converted \\\n'
                '   -t  novel_miRNA \\\n'
                '   -o  {projpath}/8.novel/novel_miRNAs \\\n'
                '   -u -s -W && \n\n'
                'awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped mature\\t"$0}}}} \' \\\n'
                '   {projpath}/8.novel/novel_miRNAs/novel_miRNA.mapmat.stat \\\n'
                '   > {projpath}/8.novel/novel_miRNAs/novel_miRNA.map.stat && \n\n'
                '{perlExec} {BioScripts}/8novel/stat_known_miRNA_pre.pl \\\n'
                '   {projpath}/8.novel/novel_miRNAs/miRBase.mrd \\\n'
                '   >{projpath}/8.novel/novel_miRNAs/novel_miRNA.mapref.stat  && \n'
                'awk \'{{if(NR==2){{print "Mapped star\\t"$0}}}}\'  \\\n'
                '    {projpath}/8.novel/novel_miRNAs/novel_miRNA.mapstar.stat  \\\n'
                '    >>{projpath}/8.novel/novel_miRNAs/novel_miRNA.map.stat && \n\n'
                'awk \'{{if(NR==2){{print "Mapped hairpin\\t"$0}}}}\' \\\n'
                '   {projpath}/8.novel/novel_miRNAs/novel_miRNA.mapref.stat \\\n' 
                '   >> {projpath}/8.novel/novel_miRNAs/novel_miRNA.map.stat && \n\n'
                'awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=int($i+0.5)}}printf("Mapped uniq sRNA\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"int($i+0.5))}}printf("\\n");}}}}\' \\\n'
                '   {projpath}/8.novel/novel_miRNAs/novel_miRNA.uc.stat \\\n'
                '   >>{projpath}/8.novel/novel_miRNAs/novel_miRNA.map.stat && \n\n'
                'awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=int($i+0.5)}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"int($i+0.5))}}printf("\\n");}}}}\' \\\n' 
                '   {projpath}/8.novel/novel_miRNAs/novel_miRNA.rc.stat \\\n'
                '   >>{projpath}/8.novel/novel_miRNAs/novel_miRNA.map.stat && \n\n'
                '{ViennaRNA2}/bin/RNAfold <  {projpath}/8.novel/novel_miRNAs/hairpin.fa \\\n'
                '   -noPS > {projpath}/8.novel/novel_miRNAs/hairpin.str && \n'
                'ln -sf {projpath}/8.novel/expression_analyses/novel.unmap.fas \\\n'
                '   {projpath}/8.novel/novel_miRNAs/novel_miRNA.unmap.fas && \n\n'
                '{perlExec} {BioScripts}/8novel/bwt12collapse.pl \\\n'
                '   {projpath}/8.novel/expression_analyses/novel_mapped.bwt.k1 \\\n'
                '   >{projpath}/8.novel/novel_miRNAs/novel_miRNA.map.fas && \n\n'
                '{perlExec} {BioScripts}/8novel/bwt_base_bias_plot.pl \\\n'
                '    -i {projpath}/8.novel/expression_analyses/novel_mapped.bwt.k1 && \n\n'
                'date +"%D %T ->Start novel miRNA stat "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 srnatoolscli = self.srnatoolscli,ViennaRNA2 = self.ViennaRNA2)
        shell_path = os.path.join(self.projpath,'8.novel','novel_miRNA_stat.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'novel_miRNA_stat', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='novel_miRNA_quantify')


    def known_TAS_and_novel_TAS_analysis(self):
        
        input_dir  = os.path.join(self.projpath,'9.TAS/input')
        make_dir(input_dir)        

        code = (
                'date +"%D %T ->Start known TAS and novel TAS analysis " && \\\n'
                'echo ==================Fetch known_TAS ================== \n'
                'cd {projpath}/9.TAS && \\\n'
                '{perlExec} {BioScripts}/8Plus_plant_TAS/known.tas_v2.pl \\\n'
                '   {refer} \\\n'
                '   {known_TAS_db} && \n'
                'echo ==================Fetch novel_TAS ================== \n'
                'export PATH={srnatoolscli}/bin:$PATH \n'
                '{perlExec} {BioScripts}/8Plus_plant_TAS/predict.tas.pl \\\n'
                '   {refer} \\\n'
                '   {projpath}/1.QC && \n'
                'echo ==================CAT known_TAS novel_TAS ================== \n'
                'cat {projpath}/9.TAS/known_TAS/known_TAS.blastn.fa \\\n'
                '    {projpath}/9.TAS/phase.out/novel_TAS.phased.fa \\\n'
                '    >{input_dir}/TAS.fna && \n'
                '{bowtie1}/bowtie-build {input_dir}/TAS.fna \\\n'
                '    {input_dir}/TAS.fna && \n'
                'date +"%D %T ->Finish known TAS and novel TAS analysis " \n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 srnatoolscli = self.srnatoolscli,refer = self.refer,known_TAS_db = self.known_TAS_db,
                 input_dir = input_dir, bowtie1 = self.bowtie1)
        shell_path = os.path.join(self.projpath,'9.TAS','known_TAS_and_novel_TAS_analysis.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'known_TAS_and_novel_TAS_analysis', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs='known_miRNA_quantify')
 

    def TAS_mapping_and_TAS_stat(self):
        
        input_dir  = os.path.join(self.projpath,'9.TAS/input')
        output_dir = os.path.join(self.projpath,'9.TAS/output')
        make_dir(output_dir)

        code = (
                'date +"%D %T ->Start TAS mapping and TAS stat " && \\\n'
                'echo ====================== Run TAS mapping ====================== \n'
                '{bowtie1}/bowtie  \\\n'
                '   -f {input_dir}/TAS.fna \\\n'
                '      {projpath}/8.novel/novel_miRNAs/novel_miRNA.unmap.fas \\\n'
                '      {output_dir}/TAS.bwt \\\n'
                '   --un {output_dir}/TAS.unmap.fas \\\n'
                '   -p 10 -v 0 -k 1 \\\n'
                '   2>{output_dir}/TAS.mapping.log && \n'
                'echo ======================= Stat TAS mapping  ====================== \n'
                'if [[ ! -z {output_dir}/TAS.bwt  ]]; then \\\n'
                '   {perlExec} {BioScripts}/8Plus_plant_TAS/bwt12collapse.pl \\\n'
                '       {output_dir}/TAS.bwt \\\n'
                '       >{output_dir}/TAS.map.fas && \n'
                '   {perlExec} {BioScripts}/8Plus_plant_TAS/genebwt12count.pl \\\n'
                '       -i {output_dir}/TAS.bwt \\\n'
                '       -r {input_dir}/TAS.fna \\\n'
                '       -t TAS -u -s \\\n'
                '       -o {output_dir} && \n'
                '   awk \'{{if(NR==1){{print "Types\\t"$0}}else if(NR==2){{print "Mapped reference\\t"$0}}}}\' \\\n'
                '       {output_dir}/TAS.mapref.stat \\\n'
                '       >{output_dir}/TAS.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped uniq sRNA\t"total);'
                'for(i=2;i<=NF;i++){{printf("\\t"$i)}}printf("\\n");}}}}\' \\\n'
                '       {output_dir}/TAS.uc.stat \\\n'
                '       >>{output_dir}/TAS.map.stat && \n'
                '   awk \'{{if(NR==2){{for(i=2;i<=NF;i++){{total+=$i}}printf("Mapped total sRNA\\t"total);'
                'for(i=2;i<=NF;i++){{printf("\t"$i)}}printf("\\n");}}}}\' \\\n' 
                '       {output_dir}/TAS.rc.stat \\\n'
                '       >>{output_dir}/TAS.map.stat \n'
                'fi && \n'
                'date +"%D %T ->Finish TAS mapping and TAS stat "\n'
        ).format(projpath = self.projpath,BioScripts = self.BioScripts,perlExec = self.perlExec,
                 bowtie1 = self.bowtie1,input_dir = input_dir, output_dir = output_dir)

        shell_path = os.path.join(self.projpath,'9.TAS','TAS_mapping_and_TAS_stat.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'TAS_mapping_and_TAS_stat', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=['known_TAS_and_novel_TAS_analysis','novel_miRNA_stat'])
        

    def sRNA_Category(self):

        mapping_stat = ','.join(
            [ 
            sampid + ":" + '{projpath}/2.map/{sampid}/{sampid}.mapping.stat'.format(
                projpath = self.projpath,sampid = sampid) \
            for sampid in self.sample_list
            ]
        )
        
        if self.org.lower() == 'refplant':
            unsense_list = (
                            '{projpath}/3.known/known_miRNAs/known_miRNA.uc.stat'
                            ':{projpath}/3.known/known_miRNAs/known_miRNA.rc.stat,'
                            '{projpath}/4.ncRNA/output/uc.stat'
                            ':{projpath}/4.ncRNA/output/rc.stat,'
                            '{projpath}/5.repeat/repeat_result/repeat.uc.stat'
                            ':{projpath}/5.repeat/repeat_result/repeat.rc.stat,'
                            '{projpath}/6.NAT/output/NAT.uc.stat'
                            ':{projpath}/6.NAT/output/NAT.rc.stat,'
                            '{projpath}/8.novel/novel_miRNAs/novel_miRNA.uc.stat'
                            ':{projpath}/8.novel/novel_miRNAs/novel_miRNA.rc.stat,'
                            '{projpath}/9.TAS/output/TAS.uc.stat'
                            ':{projpath}/9.TAS/output/TAS.rc.stat'
            ).format(projpath = self.projpath)

            sense_list   = (
                            '{projpath}/7.gene/output/uc.stat'
                            ':{projpath}/7.gene/output/rc.stat'
            ).format(projpath = self.projpath)

            previous_job = 'TAS_mapping_and_TAS_stat'
       
        else:
            if self.org.lower() == 'refanimal': 
                unsense_list = (
                                '{projpath}/3.known/known_miRNAs/known_miRNA.uc.stat'
                                ':{projpath}/3.known/known_miRNAs/known_miRNA.rc.stat,'
                                '{projpath}/4.ncRNA/output/uc.stat'
                                ':{projpath}/4.ncRNA/output/rc.stat,'
                                '{projpath}/5.repeat/repeat_result/repeat.uc.stat'
                                ':{projpath}/5.repeat/repeat_result/repeat.rc.stat,'
                                '{projpath}/8.novel/novel_miRNAs/novel_miRNA.uc.stat'
                                ':{projpath}/8.novel/novel_miRNAs/novel_miRNA.rc.stat'                                
                ).format(projpath = self.projpath)
                
                sense_list   = (
                                '{projpath}/7.gene/output/uc.stat'
                                ':{projpath}/7.gene/output/rc.stat'
                ).format(projpath = self.projpath)

                previous_job = 'novel_miRNA_stat'
                            
            else:
                unsense_list = (
                                '{projpath}/3.known/known_miRNAs/known_miRNA.uc.stat'
                                ':{projpath}/3.known/known_miRNAs/known_miRNA.rc.stat,'
                                '{projpath}/4.ncRNA/output/uc.stat'
                                ':{projpath}/4.ncRNA/output/rc.stat,'
                                '{projpath}/8.novel/novel_miRNAs/novel_miRNA.uc.stat'
                                ':{projpath}/8.novel/novel_miRNAs/novel_miRNA.rc.stat'                                
                ).format(projpath = self.projpath)
                
                previous_job = 'novel_miRNA_stat'
                
        if 'noref' not in self.org.lower():
            code = (
                    'date +"%D %T ->Start sRNA Category" && \\\n'
                    'cd {projpath}/10.Category && \\\n'
                    '{perlExec} {BioScripts}/9category/run_category.pl \\\n'
                    '   --total   {mapping_stat} \\\n'
                    '   --unsense {unsense_list} \\\n'
                    '   --sense   {sense_list} && \\\n'
                    'date +"%D %T ->Finish sRNA Category "\n'
                ).format(projpath = self.projpath, BioScripts = self.BioScripts,perlExec = self.perlExec,
                         mapping_stat = mapping_stat, sense_list = sense_list,unsense_list = unsense_list)
        else: 
            code = (
                    'date +"%D %T ->Start sRNA Category" && \\\n'
                    'cd {projpath}/10.Category && \\\n'
                    '{perlExec} {BioScripts}/9category/run_category.pl \\\n'
                    '   --total   {mapping_stat} \\\n'
                    '   --unsense {unsense_list} && \\\n'
                    'date +"%D %T ->Finish sRNA Category "\n'
                ).format(projpath = self.projpath, BioScripts = self.BioScripts,perlExec = self.perlExec,
                         unsense_list = unsense_list,mapping_stat = mapping_stat)

        shell_path = os.path.join(self.projpath,'10.Category','sRNA_Category.sh')
        write_shell(code,shell_path)
        job_name = os.path.splitext(os.path.basename(shell_path))[0]
        atom_job = AtomJob(job_name, shell_path,1,'sRNA_Category', sched=self.sched)
        self.job.add_atom_job(atom_job)
        self.job.add_order(job_name, previous_jobs=previous_job)
        self.job.add_order(job_name, after_jobs='srna_result')
