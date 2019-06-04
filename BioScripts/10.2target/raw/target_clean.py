import sys
import os
import os.path
import argparse
import glob

parser = argparse.ArgumentParser(description="qc_del_sam")
parser.add_argument('--outdir',help='filename',required=True)
parser.add_argument('--abbr',help='abbreviation of species from miRBase',required=True)
parser.add_argument('--version',help='apply for the different pipline',default="v2.3",type=str)
argv = vars(parser.parse_args())
outdir=argv['outdir'].strip()
abbr=argv['abbr']
version=argv['version']
if version=="v2.3":
    os.chdir(outdir+"/Common")
    if len(open("all_target.xls").readlines())>10:
        print "target file is OK"
        os.chdir(outdir)
        os.system("rm hsa_mature.fa_*.fasta")
        os.chdir(outdir+"/miRanda")
        os.system("rm miranda_targets_out_*")
        os.system("rm miranda_targets_out_*.fmt")
       # os.system("rm miranda_targets_*") 
        os.system("rm miranda_targets_*.pairs.example") 
        os.system("rm miranda_targets_*.pairs")
        os.system("rm core.*")
        os.chdir(outdir+"/PITA") 
        os.system("rm PITA_*_pita_results.tab") 
        os.system("rm PITA_*_pita_results_targets.tab")
        os.system("rm core.*") 
        os.chdir(outdir+"/RNAhybrid") 
        os.system("rm RNAhybrid_miRNA_target_*_pairs")
        os.system("rm core.*")
else:
    print "the version has changed, you need rewrite delete scripts part"
