import sys
import os
import os.path
import argparse
import glob

parser = argparse.ArgumentParser(description="qc_del_sam")
parser.add_argument('--outdir',help='filename',required=True)
argv = vars(parser.parse_args())
outdir=argv['outdir'].strip()

os.chdir(outdir+"/Common")
if len(open("all_target.xls").readlines())>10:
    print "target file is OK"
    os.chdir(outdir)
    os.system("rm -rf mature.fa_*.fasta")
    os.chdir(outdir+"/miRanda")
    os.system("rm -rf miranda_targets_out_*")
    os.system("rm -rf miranda_targets_out_*.fmt")
    os.system("rm -rf miranda_targets_*.pairs.example") 
    os.system("rm -rf miranda_targets_*.pairs")
    os.system("rm -rf core.*")
    os.chdir(outdir+"/PITA") 
    os.system("rm -rf PITA_*_pita_results.tab") 
    os.system("rm -rf PITA_*_pita_results_targets.tab")
    os.system("rm -rf core.*") 
    os.chdir(outdir+"/RNAhybrid") 
    os.system("rm -rf RNAhybrid_miRNA_target_*_pairs")
    os.system("rm -rf core.*")
