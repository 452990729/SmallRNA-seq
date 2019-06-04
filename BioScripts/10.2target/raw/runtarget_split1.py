#!/usr/bin/python
import os,sys
import argparse
import os.path
import glob
import ConfigParser
import misopy.cluster_utils as cluster_utils
def generate_qsub(sh, script_dir, vf):
    cluster_script = 'qsub -cwd -l vf=%dG,p=1 -o "%s" -e "%s" -V %s' % (vf, script_dir, script_dir, sh)
    return cluster_script
parser = argparse.ArgumentParser(description='split miRNA.fa file for target')
parser.add_argument('--refRNA',required=True,help='refRNA.fa file')
parser.add_argument('--type',choices=["3utr_fly","3utr_worm","3utr_human"],required=True)
parser.add_argument('--mature',required=True,help='mature fa file')
#parser.add_argument('--n',required=True,help='n seq per splited fa file')
parser.add_argument('--outdir',required=True,help='The output dir')
argv=vars(parser.parse_args())
refRNA=argv['refRNA'].strip()
type=argv['type'].strip()
mature=argv['mature'].strip()
#n=argv['n'].strip()
mature_basename=os.path.basename(mature)
mature=os.path.abspath(mature)
outdir=argv['outdir'].strip()
dir=glob.glob(os.path.dirname(outdir))[0]
os.system('mkdir -p %s/%s'%(dir,os.path.basename(outdir)))
outdir=dir+'/'+os.path.basename(outdir)
# split circ_fa file
os.system('mkdir -p %s/miRanda'%(outdir))
os.system('mkdir -p %s/PITA'%(outdir))
os.system('mkdir -p %s/RNAhybrid'%(outdir))


script_path = os.path.dirname(sys.path[0])
script_dir = os.path.dirname(script_path)
config = ConfigParser.ConfigParser()
config.read('{}/Pipeline/config.ini'.format(script_dir))
pythonExec = config.get('srnaenv','python_v276')
perlExec =  config.get('srnaenv','perl_v5182')
miranda = config.get('software','miranda')
RNAhybrid = config.get('software','RNAhybrid')

os.system('%s %s/seq_split_v1.pl %s %s/%s 20'%(perlExec,sys.path[0],mature,outdir,mature_basename))
fasta_file_list=glob.glob('%s/*.fasta'%(outdir))
refRNA=os.path.abspath(refRNA)
fa_number=len(fasta_file_list)
job_ids_target=[]
for i in range(fa_number):
	j=int(i)+1
	each_mature=outdir+'/'+mature_basename+'_'+str(j)+'.fasta'
	f1=open(outdir+"/miRanda/"'miranda_runtarget_'+str(j)+".sh",'w')
	f1.write('''
%s  %s  %s  -sc 140 -en -10 -scale 4 -strict -out %s/miRanda/miranda_targets_out_%s \n %s %s/mi_result_fmt.pl %s/miRanda/miranda_targets_out_%s  %s/miRanda/miranda_targets_out_%s.fmt \n %s %s/mi_result.pl -i %s/miRanda/miranda_targets_out_%s  -o  %s/miRanda/miranda_targets_%s  \n awk '{print substr($1,3)"\\t"$2}' %s/miRanda/miranda_targets_%s |sort -u > %s/miRanda/miranda_targets_%s.pairs
'''%(miranda,each_mature,refRNA,outdir,str(j),perlExec,sys.path[0],outdir,str(j),outdir,str(j),perlExec,sys.path[0],outdir,str(j),outdir,str(j),outdir,str(j),outdir,str(j)))
	f1.close()
	jobcmd=generate_qsub(outdir+"/miRanda/miranda_runtarget_"+str(j)+".sh",outdir+"/miRanda/",1)
	#launch_job = cluster_utils.launch_job(jobcmd, cmd_name = 'qsub')
	job_ids_target.append(cluster_utils.launch_job(jobcmd, cmd_name = 'qsub'))

for i in range(fa_number):
	j=int(i)+1	
	each_mature=outdir+'/'+mature_basename+'_'+str(j)+'.fasta'
	f2=open(outdir+"/PITA/"'pita_runtarget_'+str(j)+".sh",'w')
	f2.write('''cd  %s/PITA/ \n %s %s/pita_prediction.pl -utr %s -mir %s -prefix PITA_%s
'''%(outdir,perlExec,sys.path[0],refRNA,each_mature,str(j)))
	f2.close()
	jobcmd=generate_qsub(outdir+"/PITA/pita_runtarget_"+str(j)+".sh",outdir+"/PITA/",1)
	job_ids_target.append(cluster_utils.launch_job(jobcmd, cmd_name = 'qsub'))

for i in range(fa_number):
	j=int(i)+1
	each_mature=outdir+'/'+mature_basename+'_'+str(j)+'.fasta'
	f3=open(outdir+"/RNAhybrid/"'RNAhybrid_runtarget_'+str(j)+".sh",'w')
	f3.write('''
 %s -s %s -t %s -q %s -e -10 -p 0.05 > %s/RNAhybrid/RNAhybrid_miRNA_target_%s_pairs
'''%(RNAhybrid,type,refRNA,each_mature,outdir,str(j)))
	f3.close()
	jobcmd=generate_qsub(outdir+"/RNAhybrid/RNAhybrid_runtarget_"+str(j)+".sh",outdir+"/RNAhybrid/",1)
	#launch_job = cluster_utils.launch_job(jobcmd, cmd_name = 'qsub')
	job_ids_target.append(cluster_utils.launch_job(jobcmd, cmd_name = 'qsub'))


cluster_utils.wait_on_jobs(job_ids_target, cluster_cmd = 'qsub')

os.system('cat %s/miRanda/miranda_targets_*.pairs > %s/miRanda/miranda_targets.pairs'%(outdir,outdir))
os.system('cat %s/miRanda/miranda_targets_out_*.fmt > %s/miRanda/miranda_targets_out.fmt'%(outdir,outdir))
os.system('cat %s/PITA/PITA_*_pita_results_targets.tab > %s/PITA/PITA_pita_results_targets.tab'%(outdir,outdir))
os.system('cat %s/RNAhybrid/RNAhybrid_miRNA_target_*_pairs > %s/RNAhybrid/RNAhybrid_miRNA_target_pairs'%(outdir,outdir))

os.system('mkdir -p %s/Common' % outdir)
common_cmd ="""%s %s/get_RNAhybrid_PITA_miRanda.pl \\
        -RNAhybrid %s/RNAhybrid/RNAhybrid_miRNA_target_pairs \\
        -PITA  %s/PITA/PITA_pita_results_targets.tab   \\
        -miRanda   %s/miRanda/miranda_targets.pairs \\""" %(perlExec,sys.path[0],outdir,outdir,outdir)
open("%s/Common/common_target.sh"%outdir,"w").write(common_cmd)
os.system('cd %s/Common/ \n sh common_target.sh'%outdir)

